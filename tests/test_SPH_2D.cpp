#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <iostream>
#include <chrono>

int main(int argc, char* argv[])
{
	SPH_main domain;
	SPH_particle::main_data = &domain;

	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	domain.place_points(domain.min_x, domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	auto first_particle = domain.particle_list.front();
	write_file("particles_0.vtp", domain.particle_list);

	int iterations = 3;
	if (argc > 1) iterations = std::stoi(argv[1]);
	const double t_max = iterations*domain.dt;
	std::cout << "dt = " << domain.dt << std::endl;

	const auto start = std::chrono::high_resolution_clock::now();
	while(domain.t < t_max) {
		domain.timestep();
		if (domain.count % domain.output_intervval == 0)
		{
			write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
		}
	}
	const auto end = std::chrono::high_resolution_clock::now();
	std::cout << "runtime [s] = " << (end-start).count() / 1e9 << endl;
	write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
	
	return 0;
}
