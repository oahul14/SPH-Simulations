#pragma once

#include <list>
#include <string>
#include "SPH_2D.h"

int write_file(const string filename, const std::list<SPH_particle>& particle_list);
