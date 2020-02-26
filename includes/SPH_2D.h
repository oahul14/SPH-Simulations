#pragma once
#include <vector>
#include <utility>
#include <iterator>
#include <cmath>
#include <iostream>
#include <string>
#include <list>

using namespace std;

class SPH_main;

class SPH_particle
{

public:
	SPH_particle();
	SPH_particle(double rho, bool bound);
	double x[2], v[2];		// position and velocity
	double rho, P;					// density and pressure
	bool boundary_particle = false;
	//double rho0 = 1000;// kg/ m^3

	static SPH_main *main_data;		// link to SPH_main class so that it can be used in calc_index

	int list_num[2];				// index in neighbour finding array
	double m;
	void set_m();
	
	
	void calc_index();
	void redef_P(); //function to update the Pressure

	bool operator==(const SPH_particle& other) const;

	const SPH_particle operator+(const SPH_particle& other) const;
	const SPH_particle operator+(const SPH_particle&& other) const;
	const SPH_particle operator*(const double dt) const;
};


class SPH_main 
{
public:
	SPH_main();

	void set_values();
	void set_stencil(bool sten);
	void initialise_grid();

	void place_points(double *min, double *max, string shape = "rectangle");
    
	vector<vector<list<SPH_particle*>>> search_grid(list<SPH_particle>& particle_list);			//allocates all the points to the search grid (assumes that index has been appropriately updated)
	std::list<SPH_particle*> neighbours(const SPH_particle& part, const vector<vector<list<SPH_particle*>>> search_grid);

	SPH_particle RHS(const SPH_particle& part, const vector<vector<list<SPH_particle*>>>& search_grid);
	std::vector<SPH_particle> offsets(std::list<SPH_particle>& particle_list);
	void timestep();

	std::pair<double, double> dvdt(const SPH_particle& p, const std::list<SPH_particle*>& neighbours);

	double drhodt(const SPH_particle& p, const std::list<SPH_particle*>& neighbours);

	double W(const double r);

	double dW(const double r);

	SPH_particle smooth(const SPH_particle& part, const list<SPH_particle*>& neighbours);


	double h;								//smoothing length
	double h_fac;
	double dx = 0.2;								//particle initial spacing
	double rho0 = 1000;// kg/ m^3
	double mu = 0.001; 
	double g = 9.81;
	//to make sure equation is stiff enough
	int gamma = 7;
	//artificial speed of sound 
	double c0 = 20;
	double min_x[2], max_x[2];
	bool stencil = false;

	double t = 0;	
	double dt;
	int count = 0;
	// set the smoothing intervel in fix number
	int smoothing_interval = 10;

	//set the ouput intervel

	int output_intervval = 10;

	int max_list[2];

	list<SPH_particle> particle_list;						//list of all the particles

    //Outer 2 are the grid, inner vector is the list of pointers in each cell
    // serach grid excludes boundary cells
};
