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


#ifndef SOIL_TYPES_H
#define SOIL_TYPES_H

#include <stdio.h>
#include "soil_config.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int    ndisc;                 // must equal NDISC
    int    deepest_root_disc;     // 1..ndisc
    int    use_ch_lookup_table;   // 1 => use LUT; 0 => analytic CH
    double dt_hours;              // usually 1.0
} SoilControl;

typedef struct {
    double dz[NDISC];
    double zc[NDISC];
} SoilGeometry;

typedef struct {
    double theta_r;               // residual saturation (m3/m3)
    double theta_sat;             // saturation (m3/m3)
    double theta_fc;              // field capacity (m3/m3)
    double theta_wp;              // wilting point (m3/m3)
    double theta_aet_eq_pet;      // AET = PET at / above this \u03b8 (m3/m3)

    double K_sat_cm_per_h;        // cm/h
    double phi_sat_cm;            // cm
    double b_exp;                 // Clapp-Hornberger exponent

    double perc_limiter_0_to_1;   // 0..1 bottom drainage limiter
    double klf_m_per_h;           // lateral removal rate constant (m/h)
} SoilParameters;

// CH lookup tables over Theta=(theta-theta_r)/(theta_sat-theta_r)
typedef struct {
    int    n;
    double lnTheta_min;
    double dlnTheta;
    double inv_dlnTheta;
    double *lnpsi;                // ln(psi[m]) length n
    double *lnK;                  // ln(K[m/h]) length n
} SoilLookupTables;

typedef struct {
     double theta_in[NDISC];
     double psi_in[NDISC];   // m
     double K_in[NDISC];     // m/h
     int    ch_lut_hint_in[NDISC];
} SoilStateIn;

typedef struct {
    double theta_out[NDISC];
    int    ch_lut_hint_out[NDISC];
} SoilStateOut;

typedef struct {
    double rain_mm_per_h;         // mm/h
    double pet_mm_per_h;          // mm/h
} SoilForcing;

typedef struct {
    // Step-integrated exchanges (m)
    double AET_by_disc_m[NDISC];
    double lateral_by_disc_m[NDISC];
    double percolation_to_gw_m;
    double rain_into_soil_m;
    double rain_excess_m;

    // Internal vertical exchanges: [0..NDISC-2] interfaces i\u2192i+1; [NDISC-1] bottom perc
    double interface_vol_m[NDISC];
    double interface_rate_m_per_h[NDISC];

    int    n_sub_used;
} SoilFluxes;

typedef struct {
    double in_rain_m;         // infiltrated
    double excess_m;          // rejected
    double perc_m;
    double AET_m;
    double lateral_m;
    double delta_storage_m;   // change in storage: sum of (theta_out-theta_in)*dz
    double residual_m;        // in - (outs) - \u0394S
} TimestepSoilVolbal;

#ifdef __cplusplus
} // extern "C"
#endif

#endif // SOIL_TYPES_H
