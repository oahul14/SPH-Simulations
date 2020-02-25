#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <iostream>

int main(void)
{
	SPH_main domain;
	SPH_particle::main_data = &domain;

	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	domain.place_points(domain.min_x, domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	auto first_particle = domain.particle_list.front();
	write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);

	const double t_max = 10*domain.dt;
	std::cout << "dt = " << domain.dt << std::endl;

	while(domain.t < t_max) {
		domain.timestep();
		write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
	}
	write_file("particles_" + std::to_string(domain.t) + ".vtp", domain.particle_list);
	
	return 0;
}
