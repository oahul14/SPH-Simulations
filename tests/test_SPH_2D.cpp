#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <iostream>
#include <cassert>
#include <chrono>

// void test_timestep(int iterations = 3, SPH_main& domain) {
// 	const double t_max = iterations*domain.dt;
// 	std::cout << "dt = " << domain.dt << std::endl;
// 	write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);

// 	while(domain.t < t_max) {
// 		domain.timestep();
// 		write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
// 	}
// 	write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
// }

void test_P(SPH_main& domain) {
		//for a single particle
	auto particle = domain.particle_list.front();
	//cout<<"domain rho0 "<<domain.rho0<<"\n";
	particle.rho = 0;
	//cout<<" rho shoudl be 0: "<<particle.rho<<"\n";
	//cout<<"Particle Pressure " <<particle.P <<"\n";
	particle.redef_P();
	//cout<<"Particle Pressure after redef " <<particle.P <<"\n";
	//cout <<"Abs:" <<abs(1 - particle.P / (-57142.9))<<"\n";
	assert (abs(1 - particle.P / (-57142.9)) < 1e-6);
	//cout <<"next:" <<abs(particle.P - (-57142.9))<<"\n";
	assert (abs(particle.P - (-57142.9)) < 0.1);
	 particle.rho = 1000;
	 particle.redef_P();
	 //cout<<"Particle Pressure before 0 " <<particle.P <<"\n";
	 assert (particle.P == 0);


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
	cout<< " Wolfram "<< (-3*r)/(2*h*h)<<"\n";
	cout <<"ours "<< domain.dW(r)<<"\n";
}

int main(int argc, char* argv[])
{
	SPH_main domain;
	SPH_particle::main_data = &domain;

	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	domain.place_points(domain.min_x, domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	
	test_P(domain);
	test_W(domain);
	test_dW(domain);
	//test_timestep(domain);

	// if (argc > 1) {
	// 	test_timestep(std::stoi(argv[1]), domain);
	// }
	// else {
	// 	//test_timestep(domain)
	// }
	
	
	return 0;
}
