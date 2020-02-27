#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <iostream>
#include <cassert>
#include <chrono>

void test_timestep(SPH_main& domain, const SPH_main::timesteppers ts, const int iterations) {
	const auto t_max = iterations*domain.dt;
	std::cout << "dt = " << domain.dt << std::endl;
	write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);

	while(domain.t < t_max) {
		domain.timestep(ts);
		write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
	}
	write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
}

void test_P(SPH_main& domain) {
	//for a single particle
	auto particle = domain.particle_list.front();

	particle.rho = 0;
	particle.redef_P();
	assert (abs(1 - particle.P / (-57142.9)) < 1e-6);
	assert (abs(particle.P - (-57142.9)) < 0.1);

	particle.rho = 1000;
	particle.redef_P();
	assert (particle.P == 0);
}

void test_dvdt(SPH_main& domain) {
	auto it = domain.particle_list.begin();
	int iter=0; 
	while (iter<5)
	{
		domain.timestep();
		iter++;
	}
	
	while(it->boundary_particle) {
		it++;
	}
	auto part_1 = *it++;
	auto part_2 = *it;

	double r_ij_1 = part_1.x[0] - part_2.x[0];
	cout<<"r_ij_1: "<<r_ij_1 <<"\n";
    double r_ij_2 = part_1.x[1] - part_2.x[1];
	cout<<"r_ij_2:  "<<r_ij_2 <<"\n";
    double dist = std::sqrt(r_ij_1 * r_ij_1 + r_ij_2 * r_ij_2);
	cout<<"dist:  "<<dist <<"\n";
	double v_ij_1 = part_1.v[0] - part_2.v[0];
	cout<<"v_ij_1: "<<v_ij_1 <<"\n";
	double v_ij_2 = part_1.v[1] - part_2.v[1];
	cout<<"v_ij_2:  "<<v_ij_2 <<"\n";
	pre_calc_values pre_calculated {dist, domain.dW(dist), r_ij_1/dist, r_ij_2/dist, v_ij_1, v_ij_2};
	cout<< domain.dvdt(part_1,part_2, pre_calculated).first<<" : " << domain.dvdt(part_1,part_2, pre_calculated).first<<"\n";
	const auto c1 = -part_2.m * (part_1.P / (part_1.rho * part_1.rho) + part_2.P / (part_2.rho * part_2.rho)) * pre_calculated.dWdr;
    const auto c2 = domain.mu * part_2.m * (1 / (part_1.rho * part_1.rho) + 1 / (part_2.rho * part_2.rho)) * pre_calculated.dWdr;
    cout<<"c1:  "<<c1 <<"\n";
	cout<<"c2:  "<<c2 <<"\n";
    auto a1 = c1 * pre_calculated.e_ij_1 +  c2 * pre_calculated.v_ij_1 / pre_calculated.dist;
    auto a2 = c1 * pre_calculated.e_ij_2 +  c2 * pre_calculated.v_ij_2 / pre_calculated.dist;
	cout<<"a1:  "<<a1 <<"\n";
	cout<<"a2:  "<<a2 <<"\n";

}


void test_W(SPH_main& domain)
{
	//cout<<"H "<<domain.h<<"\n";
	//test case fof q>2
	assert (domain.W(0.65) == 0);
	//test case for 1<p<2
	//cout<<" result2 : "<< domain.W(1.5*0.26)<<"\n";
	assert((domain.W(0.39)-0.210) <1e3);
	//test case for p<1
	//cout<<"result 3: "<< domain.W(0.13)<<endl;
	assert((domain.W(0.13) - 4.8 )<0.1);
	//6.726
}

void test_dW(SPH_main& domain)
{
	double r = 0.13;
	double h = domain.h;
	cout<<"H "<< h<<"\n";
	//cout<< " Wolfram "<< (-3*r)/(2*h*h)<<"\n";
	cout<<"q: "<<r/h<<"\n";
	cout <<"dw "<< domain.dW(r/h)<<"\n";
	double a1= 3*r*(3*r-4*h)/(4*h*h*h);
	//dw =     3*r*(3*r-4*this->h)/(4*this->h*this->h*this->h);
	double factor = 10  / (7 * M_PI * h * h);
	cout<< "the components : "<< a1<<"\n";
	cout<< "the factor: "<< factor<<"\n";
	cout<< "the function: "<< a1*factor<<"\n";
}

int main(int argc, char* argv[])
{
	SPH_main domain;
	SPH_particle::main_data = &domain;

	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	domain.place_points(domain.min_x, domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	
	//test_P(domain);
	//test_W(domain);
	test_dW(domain);
	//test_dvdt(domain);

	auto ts = SPH_main::improved_euler;
	int iterations = 3;
	if (argc > 1) {
		iterations = std::stoi(argv[1]);
	}
	if (argc > 2) {
		auto arg = string(argv[2]);
		ts = (arg == "forward_euler") ? SPH_main::forward_euler : (arg == "improved_euler") ? SPH_main::improved_euler :
			throw std::invalid_argument("invalid time stepping scheme " + arg);
	}

	test_timestep(domain, ts, iterations);
	
	return 0;
}
