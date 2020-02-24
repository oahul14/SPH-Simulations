#include "SPH_2D.h"


SPH_main *SPH_particle::main_data;

void SPH_particle::calc_index(void)
{
	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0*main_data->h));
}

SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
}

void SPH_main::set_values(void)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	max_x[0] = 20.0;
	max_x[1] = 10.0;

	dx = 1;
	
	h_fac = 1.3;
	h = dx*h_fac;
}

void SPH_main::initialise_grid(void)
{
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0*h;
		max_x[i] += 2.0*h;
        //add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0*h) + 1.0);
	}

	search_grid.resize(max_list[0]);
	for (int i=0;i<max_list[0];i++)
		search_grid[i].resize(max_list[1]);
}


void SPH_main::place_points(double *min, double *max, string shape)
{
	double x[2] = { min[0], min[1] };
	SPH_particle inner_particle;
    SPH_particle outer_particle;
    inner_particle.boundary_particle = false;
    outer_particle.boundary_particle = true;
    
    if (shape == "rectangle")
    {
        while (x[1] <= max[1])
        {
            x[0] = min[0];
            while (x[0] <= max[0])
            {
                if ((x[0] < min[0] + 2.*h) || (x[0] > max[0] - 2.*h) || (x[1] < min[1] + 2.*h) || (x[1] > max[1] - 2.*h)) {
                    for (int i = 0; i < 2; i++)
                        inner_particle.x[i] = x[i];
                    inner_particle.calc_index();
                    particle_list.push_back(inner_particle);
                    cout << inner_particle.boundary_particle;
                    x[0] += dx;
                }
                else {
                    for (int i = 0; i < 2; i++)
                        outer_particle.x[i] = x[i];
                    outer_particle.calc_index();
                    particle_list.push_back(outer_particle);
                    cout << outer_particle.boundary_particle;
                    x[0] += dx;
                }
            }
            cout << endl;
            x[1] += dx;
        }
    }
    else if (shape == "shoaling")
    {
        // beach starts at half way in x1: (x1, 2h)
        double shoaling_x1_start = (max[0] - 2.*h) / 2;
        double shoaling_x1_end = max[0] - 2.*h;
        // beach ends at (max[0] - 2h, x2)
        double shoaling_x2_start = 0;
        double shoaling_x2_end = (max[1] - 2.*h) / 5;
        
        const double shoaling_slope = (shoaling_x2_end - shoaling_x2_start) / (shoaling_x1_end - shoaling_x1_start);
        cout << shoaling_slope << endl;
        const double start[2] = { shoaling_x1_start, shoaling_x2_start };
        const double end[2] = { shoaling_x1_end, shoaling_x2_end };
        while (x[1] <= max[1])
        {
            x[0] = min[0];
            while (x[0] <= max[0])
            {
                if ((x[0] < min[0] + 2.*h) || (x[0] > max[0] - 2.*h) || (x[1] < min[1] + 2.*h) || (x[1] > max[1] - 2.*h)) {
                    for (int i = 0; i < 2; i++)
                        inner_particle.x[i] = x[i];
                    inner_particle.calc_index();
                    particle_list.push_back(inner_particle);
                    cout << inner_particle.boundary_particle;
                    x[0] += dx;
                }
                else if ((x[0] >= shoaling_x1_start) && (x[0] <= shoaling_x1_end) && (x[1] >= shoaling_x2_start) && (x[1] <= shoaling_x2_end))
                {
                    for (int i = 0; i < 2; i++)
                        outer_particle.x[i] = x[i];
                    outer_particle.calc_index();
                    particle_list.push_back(outer_particle);
                    cout << inner_particle.boundary_particle;
                    x[0] += dx;
                }
                else {
                    for (int i = 0; i < 2; i++)
                        outer_particle.x[i] = x[i];
                    outer_particle.calc_index();
                    particle_list.push_back(outer_particle);
                    cout << outer_particle.boundary_particle;
                    x[0] += dx;
                }
            }
            cout << endl;
            if ((x[1] >= start[1]) && x[1] <= end[1])
            {
                shoaling_x1_start += shoaling_slope * (end[0] - end[1]);
                shoaling_x2_start += dx;
            }
            x[1] += dx;
        }
    }
}


void SPH_main::allocate_to_grid(void)				//needs to be called each time that all the particles have their positions updated
{
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}


void SPH_main::neighbour_iterate(SPH_particle *part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
	SPH_particle *other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle

	for (int i= part->list_num[0]-1;i<= part->list_num[0] + 1;i++)
		if (i>=0 && i<max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];

						if (part != other_part)					//stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
								dn[n] = part->x[n] - other_part->x[n];

							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							if (dist < 2.*h)					//only particle within 2h
							{
								//TODO: all the interactions between the particles
								
								cout << "dn: " << dn[0] << " " << dn[1] << endl;		//Should be removed from the code - simply here for you to see that it is working
							}
						}
					}
				}
}
