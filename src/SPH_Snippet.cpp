#include "C:\Users\gc2016\OneDrive - Imperial College London\ACSE\ACSE-4.3\acse-4-sph-morar\includes\SPH_2D.h"
#include "C:\Users\gc2016\OneDrive - Imperial College London\ACSE\ACSE-4.3\acse-4-sph-morar\includes\file_writer.h"
#include "C:\Users\gc2016\OneDrive - Imperial College London\ACSE\ACSE-4.3\acse-4-sph-morar\includes\navier_stokes.h"
/*
According to Microsoft:
Math Constants are not defined in Standard C/C++. 
To use them, you must first define _USE_MATH_DEFINES and then include cmath
*/
#define _USE_MATH_DEFINES
#include<cmath>
#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

SPH_main domain;

double pressure(const double rho, const double rho_0, const double c0, const double gamma, const double B)
{//implementtaion of the Tait Equition
	//relates Pressure to density
	double B = (rho_0 * pow(c0, 2) )/ gamma;
	return B * ((rho_0 / rho) - 1);
}

double W(const double r, const double h)
{
	/*
	CUBIC SPLINE: The smoothing Kernel function
	intergrates to 1 when integrated over its area
	r: distance between particles
	*/
	double q = r / h;
	double w;
	if (q < 1) { w = 1 - 1.5 * pow(q, 2) + 0.75 * pow(q,3); }
	else if (1 <= q || q < 2) { w = 0.25 * pow((2 - q), 3);}
	else { w = 0; }
	return 10*w/(7*M_PI*pow(h,2));
}

double dW(const double r, const double h)
{//differential of the cubic spline
	double q = r / h;
	double dw;
	if (q <= 1) { dw = -3 * q + (9 / 4) * pow(q, 2); }
	else { dw = -0.75 * pow((2 - q),2);}
	return 10 * dw / (7 * M_PI * pow(h, 2));
}

int main1(void)
{
	//domain.set_values();										//Set simulation parameters
	//domain.initialise_grid();									//initialise simulation grid

	//domain.place_points(domain.min_x,domain.max_x);				//places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain

	//domain.allocate_to_grid();									//needs to be called for each time step

	//domain.neighbour_iterate(&domain.particle_list[100]);		//finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle

	//write_file("example.vtp", &domain.particle_list);
	


	return 0;
}
