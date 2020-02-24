#pragma once
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

class SPH_main;

class SPH_particle
{
public:		
	double x1, x2, v1, v2;			// position and velocity
	double rho, P;					// density and pressure
	bool boundary_particle;

	static SPH_main *main_data;		// link to SPH_main class so that it can be used in calc_index

	int list_num[2];				// index in neighbour finding array
	double m = pow(main_data->dx,2)*main_data->rho0;

	void calc_index();

	bool operator==(const SPH_particle& other);
};


class SPH_main
{
public:
	SPH_main();

	void set_values(void);
	void initialise_grid(void);

	void place_points(double *min, double *max);

	void allocate_to_grid(void);			//allocates all the points to the search grid (assumes that index has been appropriately updated)

	void neighbour_iterate(SPH_particle *part);

	double h;								//smoothing length
	double h_fac;
	double dx = 0.2;								//particle initial spacing
	double rho0 = 1000;// kg/ m^3
	double mu = 0.001; 

	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell
};
