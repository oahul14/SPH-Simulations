#include "SPH_2D.h"
#include "file_writer.h"
#include <string>
#include <iostream>
#include <cassert>
#include <chrono>

void test_timeloop(SPH_main& domain, double tmax = 1, double tprec = 0.02, SPH_main::timesteppers ts = SPH_main::AB2) {
	int counter = 0;
	
	const auto start = std::chrono::high_resolution_clock::now();
	double file_write_time = 0;
	while (domain.t < tmax) {
		while(domain.t < counter * tprec) {
			if (domain.previous_offsets.empty() && ts == SPH_main::AB2) {
				domain.timestep(SPH_main::improved_euler);
			} else {
				domain.timestep(ts);
			}
		}
		counter++;
		const auto file_start = std::chrono::high_resolution_clock::now();
		write_file("particles_" + std::to_string(domain.t/100) + ".vtp", domain.particle_list);
		const auto file_end = std::chrono::high_resolution_clock::now();
		file_write_time += (file_end - file_start).count() / 1e6;
	}
	const auto end = std::chrono::high_resolution_clock::now();
	std::cout << "total iterations = " << domain.count << std::endl;
	std::cout << "avg runtime per iteration [ms] = " << ((end - start).count() / 1e6 - file_write_time) / domain.count << std::endl;
	std::cout << "avg file writing time [ms] = " << file_write_time / domain.count << std::endl;
}

int main(int argc, char* argv[]) {
    SPH_main domain;
	SPH_particle::main_data = &domain;									//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	domain.place_points(domain.min_x, domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
	
    if (argc < 3) {
        std::cout << "Usage:\nPositional arguments:\n\t<tmax>\tsimulated timespan\n\t<tprec>\tinterval for saving the state to file\nOptional arguments:\n\t<time stepping scheme>\timproved_euler (default) / forward_euler / AB2" << std::endl;
        if (argc == 2 && string(argv[1]) == "-h") return 0;
        return 1; 
    }

	SPH_main::timesteppers ts = SPH_main::forward_euler;
	const double tmax = std::stod(argv[1]);
    const double tprec = std::stod(argv[2]);
	if (argc > 3) {
		auto arg = string(argv[3]);
		ts = (arg == "forward_euler") ? SPH_main::forward_euler : (arg == "improved_euler") ? SPH_main::improved_euler : (arg == "AB2") ? SPH_main::AB2 :
			throw std::invalid_argument("invalid time stepping scheme " + arg);
	}

    test_timeloop(domain, tmax, tprec, ts);

	// test_timeloop(domain);

    return 0;
}