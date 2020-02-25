#pragma once
#include <vector>
#include <utility>
#include <iterator>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

class SPH_main;

class SPH_particle
{
public:		
	double x[2], v[2];			// position and velocity
	double rho, P;					// density and pressure
	bool boundary_particle;

	static SPH_main *main_data;		// link to SPH_main class so that it can be used in calc_index

	int list_num[2];				// index in neighbour finding array

	void calc_index();
};


class SPH_main
{
public:
	SPH_main();

	void set_values(void);
	void set_stencil(bool sten);
	void initialise_grid(void);

	void place_points(double *min, double *max, string shape);
    
	void allocate_to_grid(void);			//allocates all the points to the search grid (assumes that index has been appropriately updated)

	void neighbour_iterate(SPH_particle *part);

	bool stencil;							// stencil method or not: false by default
	double h;								//smoothing length
	double h_fac;
	double dx;								//particle initial spacing

	double min_x[2], max_x[2];
    //dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell
    // serach grid excludes boundary cells
};
