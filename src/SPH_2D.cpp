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

#include "../includes/SPH_2D.h"

#include <algorithm>
#include <list>
#include <iostream>
#include <cassert>

using namespace std;

SPH_main *SPH_particle::main_data;
double SPH_particle::B;

SPH_particle::SPH_particle()
{
	//update the Pressure evrytime you instasiate a particle
    this->rho = main_data->rho0;
    this->v[0] = this->v[1] = 0;
	this->redef_P();
}

SPH_particle::SPH_particle(double rho, bool bound)
{
    this->rho = rho;
    this->v[0] = this->v[1] = 0;
    this->set_m();
    this->redef_P();
    this->boundary_particle = bound;
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
	this->P = this->B * ((this->rho / this->main_data->rho0) - 1);
}

const offset offset::operator+(const offset& other) const {
	offset sum;
	sum.drho = this->drho + other.drho;
    sum.dx0 = this->dx0 + other.dx0;
    sum.dx1 = this->dx1 + other.dx1;
    sum.dv0 = this->dv0 + other.dv0;
    sum.dv1 = this->dv1 + other.dv1;

	return sum;
}

const offset offset::operator+(const offset&& other) const {
	offset sum = other;
	sum.drho = this->drho + other.drho;
    sum.dx0 = this->dx0 + other.dx0;
    sum.dx1 = this->dx1 + other.dx1;
    sum.dv0 = this->dv0 + other.dv0;
    sum.dv1 = this->dv1 + other.dv1;

	return sum;
}

const offset offset::operator*(const double dt) const {
	offset product;
	product.drho = this->drho * dt;
    product.dx0 = this->dx0 * dt;
    product.dx1 = this->dx1 * dt;
    product.dv0 = this->dv0 * dt;
    product.dv1 = this->dv1 * dt;

	return product;
}

void offset::operator+=(const offset& other) {
	this->drho += other.drho;
    this->dx0 += other.dx0;
    this->dx1 += other.dx1;
    this->dv0 += other.dv0;
    this->dv1 += other.dv1;
}

void SPH_particle::operator+=(const offset& delta) {
    this->rho += delta.drho;
    this->redef_P();
    if (!this->boundary_particle) {
        this->x[0] += delta.dx0;
        this->x[1] += delta.dx1;
        this->v[0] += delta.dv0;
        this->v[1] += delta.dv1;
    }
}

SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
    SPH_particle::B = (this->rho0 * pow(this->c0, 2)) / this->gamma;
}

void SPH_main::set_values(void)
{
	this->min_x[0] = 0.0;
	this->min_x[1] = 0.0;

	this->max_x[0] = 20.0;
	this->max_x[1] = 10.0;

	this->dx = 0.2;
	
	this->h_fac = 1.3;
	this->h = this->dx * this->h_fac;
    this->dt = 0.1*this->h/this->c0;
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
}


void SPH_main::place_points(double *min, double *max, string shape)
{
	double x[2] = { min[0], min[1] };

	SPH_particle water_particle(rho0, false);
    SPH_particle bound_particle(rho0, true);

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
                    particle_list.push_back(bound_particle);
                    x[0] += dx;
                }
                else if ((x[0] <= 3 && x[1] <= 5) || (3 <= x[0] && x[1] <= 2)) 
                {
                    for (int i = 0; i < 2; i++)
                        water_particle.x[i] = x[i];
                    particle_list.push_back(water_particle);
                    x[0] += dx;
                }
                else 
                {
                    x[0] += dx;
                }
            }
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
                    x[0] += dx;
                }
                else if ((x[0] >= shoaling_x1_start) && (x[0] <= shoaling_x1_end) && (x[1] >= shoaling_x2_start) && (x[1] <= shoaling_x2_end))
                {
                    for (int i = 0; i < 2; i++)
                        bound_particle.x[i] = x[i];
                    bound_particle.calc_index();
                    particle_list.push_back(bound_particle);
                    x[0] += dx;
                }
                else if ((x[0] <= 3 && x[1] <= 5) || (3 <= x[0] && x[1] <= 2)) 
                {
                    for (int i = 0; i < 2; i++)
                        water_particle.x[i] = x[i];
                    water_particle.calc_index();
                    particle_list.push_back(water_particle);
                    x[0] += dx;
                }
                else 
                {
                    x[0] += dx;
                }
            }
            // cout << endl;
            if ((x[1] >= start[1]) && x[1] <= end[1])
            {
                shoaling_x1_start += shoaling_slope * (end[0] - start[0]);
                shoaling_x2_start += dx;
            }
            x[1] += dx;
        }
    }
}

//needs to be called each time that all the particles have their positions updated
vector<vector<list<pair<SPH_particle*, list<SPH_particle*>::size_type>>>> SPH_main::search_grid(list<SPH_particle>& particle_list) {

	vector<vector<list<pair<SPH_particle*, list<SPH_particle*>::size_type>>>> search_grid(max_list[0], vector<list<pair<SPH_particle*, list<SPH_particle*>::size_type>>>(max_list[1]));

    list<SPH_particle*>::size_type pos = 0;
	for (auto& p : particle_list) {
		p.calc_index();
		search_grid[p.list_num[0]][p.list_num[1]].emplace_back(&p, pos++);
	}

	return search_grid;
}

vector<offset> SPH_main::calculate_offsets(list<SPH_particle>& particles) {

    auto search_grid = this->search_grid(particles);

    double dist;   // distance between particles
    double r_ij_1, r_ij_2; // unit vector between particles
    double v_ij_1, v_ij_2;

    vector<offset> offsets(particles.size());
    
    vector<offset>::size_type this_pos = 0;
    for (const auto& part : particles) {
        search_grid[part.list_num[0]][part.list_num[1]].pop_front();

        for (int i = max(part.list_num[0] - 1, 0); i < min(part.list_num[0] + 2, this->max_list[0]); i++) {
            for (int j = max(part.list_num[1] - 1, 0); j < min(part.list_num[1] + 2, this->max_list[1]); j++) {
                for (const auto [other_part, other_pos] : search_grid[i][j]) {
                    r_ij_1 = part.x[0] - other_part->x[0];
                    r_ij_2 = part.x[1] - other_part->x[1];
                    dist = std::sqrt(r_ij_1 * r_ij_1 + r_ij_2 * r_ij_2);
                    if (dist < 2 * this->h) {
                        v_ij_1 = part.v[0] - other_part->v[0];
                        v_ij_2 = part.v[1] - other_part->v[1];
                        this->max_vij2 = max(this->max_vij2, v_ij_1*v_ij_1 + v_ij_2*v_ij_2);

                        pre_calc_values pre_calculated {dist, this->dW(dist), r_ij_1/dist, r_ij_2/dist, v_ij_1, v_ij_2};

                        const auto [offset_part, offset_other_part] = this->calculate_offset(part, *other_part, pre_calculated);
                        offsets[other_pos] += offset_other_part;
                        offsets[this_pos] += offset_part;
                    }
                }
            }
	    }
        offsets[this_pos].dx0 = part.v[0];
        offsets[this_pos].dx1 = part.v[1];
        offsets[this_pos].dv1 -= this->g;

        this->max_ai2 = max(this->max_ai2, offsets[this_pos].dv0*offsets[this_pos].dv0 + offsets[this_pos].dv1*offsets[this_pos].dv1);

        this_pos++;
    }

    return offsets;
}

pair<offset, offset> SPH_main::calculate_offset(const SPH_particle& p_i, const SPH_particle& p_j, const pre_calc_values& vals) {
    offset offset_part, offset_neighbour;

    const auto [dv0, dv1] = this->dvdt(p_i, p_j, vals);
    offset_part.dv0 = dv0;
    offset_part.dv1 = dv1;
    offset_neighbour.dv0 = -dv0;
    offset_neighbour.dv1 = -dv1;

    offset_part.drho = offset_neighbour.drho = this->drhodt(p_i, p_j, vals);

    return make_pair(offset_part, offset_neighbour);
}


void SPH_main::timestep() 
{   
	this->max_rho = std::max_element(this->particle_list.cbegin(), this->particle_list.cend(),
                                     [](const SPH_particle& p1, const SPH_particle& p2){ return p1.rho < p2.rho; })->rho;
    this->max_vij2 = 0;
    this->max_ai2 = 0;

    const auto offsets = this->calculate_offsets(this->particle_list);

    const auto dt_cfl = this->h/sqrt(this->max_vij2);
    const auto dt_F = sqrt(this->h/sqrt(this->max_ai2));
    const auto dt_A = this->h/(this->c0*sqrt(pow(this->max_rho/this->rho0,this->gamma-1)));
    this->dt = this->Ccfl * min(dt_cfl, min(dt_F, dt_A));

	// forward euler
	auto particle_list_it = this->particle_list.begin();
	auto offsets_it = offsets.cbegin();
	while(particle_list_it != this->particle_list.end()) {
		*particle_list_it++ +=  *offsets_it++ * this->dt;
	}

    /* smoothing */
    // if (this->count > 0 && this->count % this->smoothing_interval == 0)
    // {
    //     list<SPH_particle> smoothed_state;
    //     const auto search_grid = this->search_grid(particle_list);
        
    //     for (const auto& p : this->particle_list) {
    //         smoothed_state.push_back(this->smooth(p, this->neighbours(p, search_grid)));
    //     }
    //     assert(smoothed_state.size() == this->particle_list.size());

    //     this->particle_list.swap(smoothed_state);
    // }

	// improved euler
	// const auto offsets_1 = this->offsets(this->particle_list);
	// auto next_state_star = this->particle_list;

	// auto next_state_star_it = next_state_star.begin();
	// auto offsets_1_it = offsets_1.cbegin();
	// while(next_state_star_it != next_state_star.end()) {
	// 	*next_state_star_it = *next_state_star_it + (*offsets_1_it * this->dt);
	// 	next_state_star_it->redef_P();
	// 	next_state_star_it++;
	// 	offsets_1_it++;

	// const auto offsets_2 = this->offsets(next_state_star);

	// auto particle_list_it = this->particle_list.begin();
	// offsets_1_it = offsets_1.cbegin();
	// auto offsets_2_it = offsets_2.cbegin();

	// while(particle_list_it != this->particle_list.end()) {

	// 	*particle_list_it = *particle_list_it + (*offsets_1_it + *offsets_2_it) * (0.5*this->dt);
	// 	particle_list_it->redef_P();

	// 	particle_list_it++;
	// 	offsets_1_it++;
	// 	offsets_2_it++;
	// }

    this->t += this->dt;
    this->count++;
}

double SPH_main::W(const double r)
{
	/*
	CUBIC SPLINE: The smoothing Kernel function
	intergrates to 1 when integrated over its area
	r: distance between particles
	*/
    assert(r >= 0);

	double q = r / this->h;
	double w = 0;
	if (q < 1) {
        w = 1 - 1.5 * q*q + 0.75 * q*q*q;
    }
	else if (q < 2) {
        q = 2-q;
        w = 0.25 * q*q*q;
    }

	return 10 * w / (7 * M_PI * this->h*this->h);
}

//differential of the cubic spline
double SPH_main::dW(const double r)
{
	double q = r / this->h;
	double dw = 0;
	if (q < 1) {
        dw = -3 * q + 2.25 * q*q;
    }
	else if (q < 2) {
        dw = -0.75 * (2 - q)*(2 - q);
    }
	return 10 * dw / (7 * M_PI * this->h * this->h);
}

std::pair<double, double> SPH_main::dvdt(const SPH_particle& p_i, const SPH_particle& p_j, const pre_calc_values& vals)
{
	double c1 = -p_j.m * (p_i.P / (p_i.rho * p_i.rho) + p_j.P / (p_j.rho * p_j.rho)) * vals.dWdr;
    double c2 = this->mu * p_j.m * (1 / (p_i.rho * p_i.rho) + 1 / (p_j.rho * p_j.rho)) * vals.dWdr;

    return make_pair(c1 * vals.e_ij_1 +  c2 * vals.v_ij_1 / vals.dist,
                     c1 * vals.e_ij_2 +  c2 * vals.v_ij_2 / vals.dist);
}

double SPH_main::drhodt(const SPH_particle& p_i, const SPH_particle& p_j, const pre_calc_values& vals)
{
	return p_j.m * vals.dWdr * (vals.v_ij_1 * vals.e_ij_1 + vals.v_ij_2 * vals.e_ij_2);
}

SPH_particle SPH_main::smooth(const SPH_particle& part, const list<pair<SPH_particle*, pre_calc_values>>& neighbours)
{
    auto smoothed = part;
    double w;
    double sum_w = 0;
    double sum_wdrho = 0;

    for (const auto& [part, vals] : neighbours)
    {
        w = this->W(vals.dist);
        sum_w += w;
        sum_wdrho += w / part->rho;
    }
    smoothed.rho = sum_w / sum_wdrho;

    return smoothed;
}