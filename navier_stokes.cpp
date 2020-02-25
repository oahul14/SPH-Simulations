#include "navier_stokes.h"
#include <cmath>
#include "SPH_2D.h"

double dW(const double r, const double h)
{
    return 2.2;
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