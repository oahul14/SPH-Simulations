#include <utility>
#include <list>
#include "SPH_2D.h"


std::pair<double, double> dvdt(const SPH_particle& p, const std::vector<SPH_particle>& neighbours);

double drhodt(const SPH_particle& p, const std::vector<SPH_particle>& neighbours);

double W(const double r, const double h);

double dW(const double r, const double h);