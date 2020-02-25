#include "navier_stokes.h"
#include <cmath>
#include "SPH_2D.h"
/*
According to Microsoft:
Math Constants are not defined in Standard C/C++.
To use them, you must first define _USE_MATH_DEFINES and then include cmath
in some libraries the M_PI is not include so we included the #ifndef
*/
#define _USE_MATH_DEFINES
#include<cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double pressure()
{
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
    if (q < 1) { w = 1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3); }
    else if (1 <= q || q < 2) { w = 0.25 * pow((2 - q), 3); }
    else { w = 0; }
    return 10 * w / (7 * M_PI * pow(h, 2));
}

double dW(const double r, const double h)
{//differential of the cubic spline
    double q = r / h;
    double dw;
    if (q <= 1) { dw = -3 * q + (9 / 4) * pow(q, 2); }
    else { dw = -0.75 * pow((2 - q), 2); }
    return 10 * dw / (7 * M_PI * pow(h, 2));
}

std::pair<double, double> dvdt(const SPH_particle& p, const std::vector<SPH_particle>& neighbours)
{
    std::pair<double,double> a(0.0,0.0);
    for (const auto& i : neighbours)
    {
        if(i==p) continue;
        else
        {
            double r_ij_1 = p.x[0] - i.x[0];
            double r_ij_2 = p.x[1] - i.x[1];
            double dist = std::sqrt(std::pow(r_ij_1,2)+ std::pow(r_ij_2,2));
            double v_ij_1 = p.v[0] - i.v[0];
            double v_ij_2 = p.v[1] - i.v[1];
            double e_ij_1 = r_ij_1/dist;
            double e_ij_2 = r_ij_2/dist;

            double dwdr = dW(dist, p.main_data->h);
            
            // x direction
            double a1 = -i.m*((p.P/pow(p.rho,2)) + (i.P/pow(i.rho,2)))*dwdr*e_ij_1 + p.main_data->mu*(i.m*(1/pow(p.rho,2)+1/pow(i.rho,2))*dwdr*(v_ij_1/dist));

            // y direction
            double a2 = -i.m*((p.P/pow(p.rho,2)) + (i.P/pow(i.rho,2)))*dwdr*e_ij_2 + p.main_data->mu*(i.m*(1/pow(p.rho,2)+1/pow(i.rho,2))*dwdr*(v_ij_2/dist));

            
            a.first = a.first + a1;
            a.second = a.second + a2;

        }
    }
    return a; 
}
double drhodt(const SPH_particle& p, const std::list<SPH_particle>& neighbours)
{
    double D;
    for (const auto& i : neighbours)
    {
        double r_ij_1 = p.x[0] - i.x[0];
        double r_ij_2 = p.x[1] - i.x[1];
        double dist = std::sqrt(std::pow(r_ij_1,2)+ std::pow(r_ij_2,2));
        double v_ij_1 = p.v[0] - i.v[0];
        double v_ij_2 = p.v[1] - i.v[1];
        double e_ij_1 = r_ij_1/dist;
        double e_ij_2 = r_ij_2/dist;

        double dwdr = dW(dist, p.main_data->h);

        D = D + i.m*dwdr*(v_ij_1*e_ij_1 + v_ij_2*e_ij_2);
    }
    return D;

}