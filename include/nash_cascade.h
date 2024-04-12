#ifndef _NASH_H
#define _NASH_H

#define MAX_NUM_NASH_CASCADE  3

// this data structure describes runoff using Nash Cascade model in the subsurface and on the surface
struct nash_cascade_parameters {
  int    N_nash;           // Number of Nash cascade reservoirs; [-]
  double K_nash;           // Fraction of storage per hour that moves from one reservoir to next (time constant); [1/hour]
  int    nsubsteps;        // the number of substeps that each dt is divided into
  double *nash_storage;    // storage array nash cascade reservoirs; [m]
};

extern double nash_cascade(double flux_lat_m,int num_lateral_flow_nash_reservoirs,
                           double K_nash,double *nash_storage_arr);

double nash_cascade_surface_runoff(double runoff_m, struct nash_cascade_parameters *nash_params);

#endif
