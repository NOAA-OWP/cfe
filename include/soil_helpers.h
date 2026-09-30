/*
 * NOAA-OWP/cfe - Version 3 of Conceptual Functional Equivalent to the stormflow/runoff
 *                generation components of the NOAA/NWS National Water Model version 3.1 
 *                and earlier
 *
 * Originally conceived and developed by: 
 *         Fred L. Ogden, Chief Scientist, NOAA/NWS 
 *         Office of Water Prediction, Tuscaloosa, AL
 *
 */


#ifndef SOIL_HELPERS_H
#define SOIL_HELPERS_H

#include <stddef.h>
#include "cfe_soil_discrete.h"
#include "cfe_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------- CH analytic relations (psi in m, K in m/h) ---------- */
double psi_from_theta(double theta, double theta_sat,
                      double phi_sat_cm, double b_exp);

double theta_from_psi(double psi_m, double theta_sat,
                      double phi_sat_cm, double b_exp);

double K_from_theta(double theta, double theta_sat,
                    double K_sat_cm_per_h, double b_exp);

/* ---------- LUT build / free ---------- */
int  soil_build_ch_lut(int n_points, double Theta_min,
                       double theta_r, double theta_sat, double b_exp,
                       double phi_sat_cm, double K_sat_cm_per_h,
                       SoilLookupTables *lut_out);

void soil_free_ch_lut(SoilLookupTables *lut);

/* Evaluate psi,K from theta using LUT if provided; updates hint in place */
void eval_theta_to_props(const SoilLookupTables *lut,
                         double theta,
                         double theta_r, double theta_sat,
                         double phi_sat_cm, double K_sat_cm_per_h, double b_exp,
                         int *hint, double *psi_m_out, double *K_m_per_h_out);

/* Vector wrapper (updates hint array in place) */
void compute_props_with_option_stateless(const SoilLookupTables *lut,
                                         const double theta[NDISC],
                                         double psi_m_out[NDISC],
                                         double K_m_per_h_out[NDISC],
                                         double theta_r, double theta_sat,
                                         double K_sat_cm_per_h, double phi_sat_cm, double b_exp,
                                         int hint_inout[NDISC]);

/* ---------- Fluxes ---------- */
double flux_DB_pair(double psi_up, double K_up,
                    double psi_dn, double K_dn,
                    double dz_up,  double dz_dn);      // m/h, downward positive

/* Choose substeps using dt_hours (generic across NDISC) */
int choose_n_sub_dt(double rain_mm_per_h,
                    double theta1, double theta_sat, double dz1,
                    double q12_0, double q23_0, double q34_0, // legacy 4-disc; keep for tests
                    double dzmin, double dt_hours);

/* Lateral removal (m removed this substep) */
double remove_lateral_to_subsurface_nash_substep(double theta[NDISC], const double dz[NDISC],
                                                 double theta_fc, double theta_sat,
                                                 double rate_const_m_per_h, double dt_sub,
                                                 double removed_by_disc_accum_m[NDISC]);

/* ---------- Hydrostatic initialization ---------- */
void initialize_hydrostatic_from_storage(double soil_depth_m,
                                         double theta_sat, double phi_sat_cm, double b_exp,
                                         const double zc[NDISC],
                                         double target_storage_m,
                                         double *zwt_out,
                                         double theta_out[NDISC]);



#ifdef __cplusplus
} // extern "C"
#endif

#endif // SOIL_HELPERS_H
