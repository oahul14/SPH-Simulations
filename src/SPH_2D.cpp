#include "../includes/SPH_2D.h"
/*
According to Microsoft:
Math Constants are not defined in Standard C/C++.
To use them, you must first define _USE_MATH_DEFINES and then include cmath
in some libraries the M_PI is not include so we included the #ifndef
*/
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


SPH_main *SPH_particle::main_data;

SPH_particle::SPH_particle()
{
	//update the Pressure evrytime you instasiate a particle
    rho = main_data->rho0;
	redef_P();
}

SPH_particle::SPH_particle(double rho, bool bound)
{
    rho = rho;
    set_m();
    redef_P();
    boundary_particle = bound;
}

void SPH_particle::calc_index(void)
{
	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0*main_data->h));
}

void SPH_particle::set_m()
{
	m = pow(main_data->dx, 2) * (main_data->rho0);
}
void SPH_particle::redef_P()
{
	//implementtaion of the Tait Equition
	//relates Pressure to density
	double B = (main_data->rho0 * pow(main_data->c0, 2)) / main_data->gamma;
	P = B * ((main_data->rho0 / rho) - 1);
}

SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
    stencil = false;
}

void SPH_main::set_values(void)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	max_x[0] = 20.0;
	max_x[1] = 10.0;

	dx = 0.2;
	
	h_fac = 1.3;
	h = dx*h_fac;
}

void SPH_main::set_stencil(bool sten) { stencil = sten; }

void SPH_main::initialise_grid(void)
{
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0*h;
		max_x[i] += 2.0*h;
        //add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0*h) + 1.0);
	}

	search_grid.resize(max_list[0]);
	for (int i=0;i<max_list[0];i++)
		search_grid[i].resize(max_list[1]);
}


void SPH_main::place_points(double *min, double *max, string shape)
{
    cout << min[0] << " " << min[1] << endl;
	double x[2] = { min[0], min[1] };
	SPH_particle water_particle(rho0, false);
    SPH_particle bound_particle(rho0, true);
    SPH_particle empty_particle(0, false);
    
    if (shape == "rectangle")
    {
        while (x[1] <= max[1])
        {
            x[0] = min[0];
            while (x[0] <= max[0])
            {
                if ((x[0] < min[0] + 2.*h) || (x[0] > max[0] - 2.*h) || (x[1] < min[1] + 2.*h) || (x[1] > max[1] - 2.*h)) 
                {
                    for (int i = 0; i < 2; i++)
                        bound_particle.x[i] = x[i];
                    bound_particle.calc_index();
                    particle_list.push_back(bound_particle);
                    // cout << particle_list[particle_list.size() - 1].boundary_particle;
                    cout << "B ";
                    x[0] += dx;
                }
                else if ((x[0] <= 3 && x[1] <= 5) || (3 <= x[0] && x[1] <= 2)) 
                {
                    for (int i = 0; i < 2; i++)
                        water_particle.x[i] = x[i];
                    water_particle.calc_index();
                    particle_list.push_back(water_particle);
                    // cout << particle_list[particle_list.size() - 1].boundary_particle;
                    cout << "W ";
                    x[0] += dx;
                }
                else 
                {
                    for (int i = 0; i < 2; i++)
                        empty_particle.x[i] = x[i];
                    empty_particle.calc_index();
                    particle_list.push_back(empty_particle);
                    // cout << particle_list[particle_list.size() - 1].boundary_particle;
                    cout << "E ";
                    x[0] += dx;
                }
            }
            cout << endl;
            x[1] += dx;
        }
    }
    else if (shape == "shoaling")
    {
        double beach_height_ratio = 0.2;
        double beach_length_ratio = 0.5;
        
        // beach starts at half way in x1: (x1, 2h)
        double shoaling_x1_start = (max[0] - 2.*h) - beach_length_ratio * (max[0] - 2.*h);
        double shoaling_x1_end = max[0] - 2.*h;
        // beach ends at (max[0] - 2h, x2)
        double shoaling_x2_start = 0;
        double shoaling_x2_end = beach_height_ratio * (max[1] - 2.*h) - dx;
        
        const double shoaling_slope = dx / (shoaling_x2_end - shoaling_x2_start);
        const double start[2] = { shoaling_x1_start, shoaling_x2_start };
        const double end[2] = { shoaling_x1_end, shoaling_x2_end };
        while (x[1] <= max[1])
        {
            x[0] = min[0];
            while (x[0] <= max[0])
            {
                if ((x[0] < min[0] + 2.*h) || (x[0] > max[0] - 2.*h) || (x[1] < min[1] + 2.*h) || (x[1] > max[1] - 2.*h)) 
                {
                    for (int i = 0; i < 2; i++)
                        bound_particle.x[i] = x[i];
                    bound_particle.calc_index();
                    particle_list.push_back(bound_particle);
                    // cout <<  particle_list[particle_list.size() - 1].boundary_particle;
                    cout << "B ";
                    x[0] += dx;
                }
                else if ((x[0] >= shoaling_x1_start) && (x[0] <= shoaling_x1_end) && (x[1] >= shoaling_x2_start) && (x[1] <= shoaling_x2_end))
                {
                    for (int i = 0; i < 2; i++)
                        bound_particle.x[i] = x[i];
                    bound_particle.calc_index();
                    particle_list.push_back(bound_particle);
                    // cout <<  particle_list[particle_list.size() - 1].boundary_particle;
                    cout << "B ";
                    x[0] += dx;
                }
                else if ((x[0] <= 3 && x[1] <= 5) || (3 <= x[0] && x[1] <= 2)) 
                {
                    for (int i = 0; i < 2; i++)
                        water_particle.x[i] = x[i];
                    water_particle.calc_index();
                    particle_list.push_back(water_particle);
                    // cout << particle_list[particle_list.size() - 1].boundary_particle;
                    cout << "W ";
                    x[0] += dx;
                }
                else 
                {
                    for (int i = 0; i < 2; i++)
                        empty_particle.x[i] = x[i];
                    empty_particle.calc_index();
                    particle_list.push_back(empty_particle);
                    // cout << particle_list[particle_list.size() - 1].boundary_particle;
                    cout << "E ";
                    x[0] += dx;
                }
            }
            cout << endl;
            if ((x[1] >= start[1]) && x[1] <= end[1])
            {
                shoaling_x1_start += shoaling_slope * (end[0] - start[0]);
                shoaling_x2_start += dx;
            }
            x[1] += dx;
        }
    }
}


void SPH_main::allocate_to_grid(void)				//needs to be called each time that all the particles have their positions updated
{
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}


void SPH_main::neighbour_iterate(SPH_particle *part)					
//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
    double dist;			//distance between particles
    double dn[2];			//vector from 1st to 2nd particle

    if (!stencil) {
        SPH_particle *other_part;
        for (int i= part->list_num[0]-1;i<= part->list_num[0] + 1;i++)
            if (i>=0 && i<max_list[0])
                for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
                    if (j >= 0 && j < max_list[1])
                    {
                        for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                        {
                            other_part = search_grid[i][j][cnt];
                            if (part != other_part)					//stops particle interacting with itself
                            {
                                //Calculates the distance between potential neighbours
                                for (int n = 0; n < 2; n++)
                                    dn[n] = part->x[n] - other_part->x[n];

                                dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

                                if (dist < 2.*h)					//only particle within 2h
                                {
                                    //TODO: all the interactions between the particles
                                    //Should be removed from the code - simply here for you to see that it is working
                                    cout << "dn: " << dn[0] << " " << dn[1] << endl;	
                                    neighbours.push_back(other_part);	
                                }
                            }
                        }
                    }
    }
    
    
    // stencil
    else {
        int box_i = part->list_num[0];
        int box_j = part->list_num[1];
        vector<pair<int, int>> stencils;

        stencils.push_back(make_pair(box_i, box_j));
        stencils.push_back(make_pair(box_i + 1, box_j));
        stencils.push_back(make_pair(box_i + 1, box_j + 1));
        stencils.push_back(make_pair(box_i, box_j + 1));
        stencils.push_back(make_pair(box_i - 1, box_j + 1));
        // loop over boxes in stencils
        for (auto box : stencils)
            // make sure the neighbour box within grid
            if ((0 <= box.first < max_list[0]) && (0 <= box.second < max_list[1]))
                for (auto opart : search_grid[box.first][box.second])
                    if (part != opart)
                    {
                        for (int n = 0; n < 2; n++)
                            dn[n] = part->x[n] - opart->x[n];
                        dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
                        if (dist < 2.*h)
                        {
                            cout << "dn: " << dn[0] << " " << dn[1] << endl;
                            neighbours.push_back(opart);
                        }
                    }
    }
}

bool SPH_particle::operator==(const SPH_particle& other) const
{
	if(this->x[0] != other.x[0]) return false;
	if(this->x[1] != other.x[1]) return false;
	if(this->v[0] != other.v[0]) return false;
	if(this->v[1] != other.v[1]) return false;
	if(this->rho != other.rho) return false;
	else return true;
}

double SPH_main::W(const double r)
{
	/*
	CUBIC SPLINE: The smoothing Kernel function
	intergrates to 1 when integrated over its area
	r: distance between particles
	*/
	double q = r / h;
	double w;
	if (q < 1) { w = 1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3); }
	else if (1 <= q || q < 2) { w = 0.25 * pow((2 - q), 3); }
	else { w = 0; }
	return 10 * w / (7 * M_PI * pow(h, 2));
}

double SPH_main::dW(const double r)
{//differential of the cubic spline
	double q = r / h;
	double dw;
	if (q <= 1) { dw = -3 * q + (9 / 4) * pow(q, 2); }
	else { dw = -0.75 * pow((2 - q), 2); }
	return 10 * dw / (7 * M_PI * pow(h, 2));
}

std::pair<double, double> SPH_main::dvdt(const SPH_particle& p, const std::vector<SPH_particle>& neighbours)
{
	std::pair<double, double> a(0.0, 0.0);
	for (const auto& i : neighbours)
	{
		if (i == p) continue;
		else
		{
			double r_ij_1 = p.x[0] - i.x[0];
			double r_ij_2 = p.x[1] - i.x[1];
			double dist = std::sqrt(std::pow(r_ij_1, 2) + std::pow(r_ij_2, 2));
			double v_ij_1 = p.v[0] - i.v[0];
			double v_ij_2 = p.v[1] - i.v[1];
			double e_ij_1 = r_ij_1 / dist;
			double e_ij_2 = r_ij_2 / dist;

			double dwdr = dW(dist);

			// x direction
			double a1 = -i.m * ((p.P / pow(p.rho, 2)) + (i.P / pow(i.rho, 2))) * dwdr * e_ij_1 + p.main_data->mu * (i.m * (1 / pow(p.rho, 2) + 1 / pow(i.rho, 2)) * dwdr * (v_ij_1 / dist));

			// y direction
			double a2 = -i.m * ((p.P / pow(p.rho, 2)) + (i.P / pow(i.rho, 2))) * dwdr * e_ij_2 + p.main_data->mu * (i.m * (1 / pow(p.rho, 2) + 1 / pow(i.rho, 2)) * dwdr * (v_ij_2 / dist));


			a.first = a.first + a1;
			a.second = a.second + a2;

		}
	}
	return a;
}
double SPH_main::drhodt(const SPH_particle& p, const std::vector<SPH_particle>& neighbours)
{
	double D;
	for (const auto& i : neighbours)
	{
		double r_ij_1 = p.x[0] - i.x[0];
		double r_ij_2 = p.x[1] - i.x[1];
		double dist = std::sqrt(std::pow(r_ij_1, 2) + std::pow(r_ij_2, 2));
		double v_ij_1 = p.v[0] - i.v[0];
		double v_ij_2 = p.v[1] - i.v[1];
		double e_ij_1 = r_ij_1 / dist;
		double e_ij_2 = r_ij_2 / dist;

		double dwdr = dW(dist);

		D = D + i.m * dwdr * (v_ij_1 * e_ij_1 + v_ij_2 * e_ij_2);
	}
	return D;
}

void SPH_main::smooth(SPH_particle *part)
{
    double w = W(0.);       // calculate kernal for part itself
    double sum_w = w;
    double sum_wdrho = w / part->rho;
    for (auto opart : neighbours)
    {
        double r_ij = sqrt(pow(opart->x[0] - part->x[0], 2) + pow(opart->x[1] - part->x[1], 2));
        w = W(r_ij);
        sum_w += w;
        sum_wdrho += w / opart->rho;
    }
    part->rho = sum_w / sum_wdrho;
}