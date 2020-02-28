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

struct pre_calc_values {
	double dist, dWdr, e_ij_1, e_ij_2, v_ij_1, v_ij_2; 
};

struct offset {
	double dx0 = 0;
	double dx1 = 0;
	double dv0 = 0;
	double dv1 = 0;
	double drho = 0;

	const offset operator+(const offset& other) const;
	const offset operator+(const offset&& other) const;
	const offset operator*(const double dt) const;

	void operator+=(const offset& other);
};

class SPH_particle
{

public:
	SPH_particle();
	SPH_particle(double rho, bool bound);
	double x[2], v[2];		// position and velocity
	double rho, P;					// density and pressure
	bool boundary_particle;

	static SPH_main *main_data;		// link to SPH_main class so that it can be used in calc_index
	static double B;

	int list_num[2];				// index in neighbour finding array
	double m;
	void set_m();
	
	
	void calc_index();
	void redef_P(); //function to update the Pressure

	void operator+=(const offset& delta);
};

class SPH_main 
{
public:
	enum timesteppers {
		forward_euler, improved_euler, AB2
	};

	SPH_main();
	void initialise_grid();

	void place_points(double *min, double *max, string shape = "rectangle");
    
	//allocates all the points to the search grid (assumes that index has been appropriately updated)
	vector<vector<list<pair<SPH_particle*, list<SPH_particle*>::size_type>>>> search_grid(list<SPH_particle>& particle_list);

	vector<offset> calculate_offsets(list<SPH_particle>& particles);

	pair<offset, offset> calc_offset(const SPH_particle& p_i, const SPH_particle& p_j, const pre_calc_values& vals);
	list<offset> offsets(list<SPH_particle>& particle_list);
	void timestep(const timesteppers& ts = forward_euler);

	pair<double, double> dvdt(const SPH_particle& p_i, const SPH_particle& p_j, const pre_calc_values& vals);
	double drhodt(const SPH_particle& p_i, const SPH_particle& p_j, const pre_calc_values& vals);
	double W(const double r);
	double dW(const double r);

	void smooth(list<SPH_particle>& particles);

	void drag_back(SPH_particle& part);

	double h_fac = 1.3;
	double h;
	double dx = 0.1;								//particle initial spacing
	double rho0 = 1000;// kg/ m^3
	double mu = 0.001; 
	double g = 9.81;
	//to make sure equation is stiff enough
	int gamma = 7;
	//artificial speed of sound 
	double c0 = 20;
	// this parameter is in [0.1,0,3]
	double Ccfl = 0.2;
	double min_x[2] {0, 0};
	double max_x[2] {20, 10};

	double max_vij2, max_ai2, max_rho;

	double t = 0;	
	double dt;
	double prev_dt = 0;
	int count = 0;
	// set the smoothing intervel in fix number
	int smoothing_interval = 10;

	//set the ouput intervel

	int output_intervval = 10;

	int max_list[2];

	list<SPH_particle> particle_list;						//list of all the particles
	vector<offset> previous_offsets;
    //Outer 2 are the grid, inner vector is the list of pointers in each cell
    // serach grid excludes boundary cells
};
