#ifndef _NASH_H
#define _NASH_H

#define MAX_NUM_NASH_CASCADE    3

extern double nash_cascade(double flux_lat_m,int num_lateral_flow_nash_reservoirs,
                           double K_nash,double *nash_storage_arr);
#endif
