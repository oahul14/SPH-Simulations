#include <utility>
#include <list>

#include "SPH_2D.h"

std::pair<double, double> dvdt(const SPH_particle& p, const std::list<SPH_particle>& neighbours);

double drhodt(const SPH_particle& p, const std::list<SPH_particle>& neighbours);

double pressure(const double rho, const double rho_0, const double c0, const double gamma, const double B);

double W(const double r, const double h);

double dW(const double r, const double h);