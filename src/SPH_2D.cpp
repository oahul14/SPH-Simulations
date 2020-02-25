#include "C:\Users\gc2016\OneDrive - Imperial College London\ACSE\ACSE-4.3\acse-4-sph-morar\includes\SPH_2D.h"
/*
According to Microsoft:
Math Constants are not defined in Standard C/C++.
To use them, you must first define _USE_MATH_DEFINES and then include cmath
in some libraries the M_PI is not include so we included the #ifndef
*/
#define _USE_MATH_DEFINES
#include<cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

SPH_main *SPH_particle::main_data;

SPH_particle::SPH_particle()
{
	//update the Pressure evrytime you instasiate a particle
	redef_P();
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
}

void SPH_main::set_values(void)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	max_x[0] = 1.0;
	max_x[1] = 1.0;

	dx = 0.02;
	
	h_fac = 1.3;
	h = dx*h_fac;
}

void SPH_main::initialise_grid(void)
{
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0*h;
		max_x[i] += 2.0*h;												//add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0*h) + 1.0);
	}

	search_grid.resize(max_list[0]);
	for (int i=0;i<max_list[0];i++)
		search_grid[i].resize(max_list[1]);
}


void SPH_main::place_points(double *min, double *max)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;

	while (x[0] <= max[0])
	{
		x[1] = min[1];
		while (x[1] <= max[1])
		{
			for (int i = 0; i < 2; i++)
				particle.x[i] = x[i];

			particle.calc_index();

			particle_list.push_back(particle);

			x[1] += dx;
		}
		x[0] += dx;
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


void SPH_main::neighbour_iterate(SPH_particle *part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
	SPH_particle *other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle

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
								
								cout << "dn: " << dn[0] << " " << dn[1] << endl;		//Should be removed from the code - simply here for you to see that it is working
							}
						}
					}
				}
}
bool SPH_particle::operator==( SPH_particle& other)
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