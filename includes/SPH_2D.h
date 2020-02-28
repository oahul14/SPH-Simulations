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

	/**
	 * @brief velocity and position stored in two directions
	 */
	double x[2], v[2];		

	/**
	 * @brief density rho and pressure p
	 */
	double rho, P;					
	
	/**
	 * @brief boolean to indicate boundary particle or not 
	 */
	bool boundary_particle;

	/**
	 * @brief link to SPH_main class so that it can be used in calc_index
	 */
	static SPH_main *main_data;		


	static double B;

	/**
	 * @brief index in neighbour finding array
	 * 
	 */
	int list_num[2];	

	/**
	 * @brief mass for particle
	 */
	double m;

	/**
	 * @brief Set the m object
	 */
	void set_m();
	
	/**
	 * @brief calculate the index for grids into list_num
	 */
	void calc_index();

	/**
	 * @brief update pressure using rho: tait equation
	 */
	void redef_P(); //function to update the Pressure

	/**
	 * @brief operator overload for struct offset to update current stage
	 * 
	 * @param delta offset variable storing dxdt, dvdt and drhodt
	 */
	void operator+=(const offset& delta);
};

class SPH_main 
{
public:
	enum timesteppers {
		forward_euler, improved_euler, AB2
	};
	
	/**
	 * @brief create shapes enum for placing points: rectangle and shoaling
	 */
	enum shapes {
		rectangle, shoaling
	};

	SPH_main();
	/**s
	 * @brief initialising min_x and max_x and max number of boxes of 2h
	 */
	void initialise_grid();

	/**
	 * @brief placing water cells within inner domain; placing boundary cells for outer domain
	 * 
	 * @param min domain's minimum x and y coordinates for the entire domain (including boundary)
	 * @param max domain's maximum x and y coordinates for the entire domain (including boundary)
	 * @param shape boundary shape: including normal rectangle and shoaling
	 */
	void place_points(double *min, double *max, const shapes& shape = rectangle);
    
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

	/**
	 * @brief characteristic smoothing length
	 */
	double h;

	/**
	 * @brief the factor to determine h: set to 1.3
	 */
	double h_fac = 1.3;

	/**
	 * @brief mesh resolution
	 */
	double dx = 0.1;	

	/**
	 * @brief initial density for both water and boundary cells, in kg/m3
	 */
	double rho0 = 1000;

	/**
	 * @brief viscosity constant in navier stokes equation
	 */
	double mu = 0.001; 

	/**
	 * @brief gravitational acceleration constant
	 */
	double g = 9.81;

	/**
	 * @brief exponential for the ratio of current rho and initial rho 
	 * to provide stiff enough (slightly compressible) state
	 */
	int gamma = 7;
	
	/**
	 * @brief artificial velocity of sound
	 * selected to be larger than the speed sustained in the system
	 */
	double c0 = 20;

	/**
	 * @brief Courant–Friedrichs–Lewy condition: should be in range [0.1, 0.3], set as 0.2
	 */
	double Ccfl = 0.2;

	/**
	 * @brief minimum x and y coordinates: initialised as the size of the inner domain
	 */
	double min_x[2] {0, 0};

	/**
	 * @brief maximum x and y coordinates: initialised as the size of the inner domain
	 * 
	 */
	double max_x[2] {20, 10};

	double max_vij2, max_ai2, max_rho;

	/**
	 * @brief current time state for the domain 
	 */
	double t = 0;	

	/**
	 * @brief dynamic time step
	 * 
	 */
	double dt;

	/**
	 * @brief get track of dt at previous stage
	 */
	double prev_dt = 0;

	int count = 0;

	/**
	 * @brief set the smoothing intervel in fix number
	 */
	int smoothing_interval = 10;

	/**
	 * @brief set the ouput intervel
	 */
	int output_intervval = 10;

	int max_list[2];

	/**
	 * @brief list to store all particles
	 */
	list<SPH_particle> particle_list;

	/**
	 * @brief vector to keep track of the previous stage offsets
	 */
	vector<offset> previous_offsets;
};
