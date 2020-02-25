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
	SPH_particle();
	double x[2], v[2];		// position and velocity
	double rho, P;					// density and pressure
	bool boundary_particle;
	//double rho0 = 1000;// kg/ m^3

	static SPH_main *main_data;		// link to SPH_main class so that it can be used in calc_index

	int list_num[2];				// index in neighbour finding array
	double m;
	void set_m();
	
	
	void calc_index();
	void redef_P(); //function to update the Pressure

	bool operator==(const SPH_particle& other) const;
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

	std::pair<double, double> dvdt(const SPH_particle& p, const std::vector<SPH_particle>& neighbours);

	double drhodt(const SPH_particle& p, const std::vector<SPH_particle>& neighbours);

	double W(const double r);

	double dW(const double r);


	double h;								//smoothing length
	double h_fac;
	double dx = 0.2;								//particle initial spacing
	double rho0 = 1000;// kg/ m^3
	double mu = 0.001; 
	//to make sure equation is stiff enough
	int gamma = 7;
	//artificial speed of sound 
	double c0 = 20;
	double min_x[2], max_x[2];
	bool stencil;	

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell
    // serach grid excludes boundary cells
};
