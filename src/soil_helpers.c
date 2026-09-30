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
 

// soil_helpers.c  
// ASCII only; banners; flux positive downward; units noted on each function.
// This file contains: constants, Look Up Table (LUT) structs/build/eval, 
// Clapp-Hornberger (CH) analytic relations, hydrostatic initializers, 
// Darcy-Buckingham (DB) flux calculation, code to calculate the number of
// sub timesteps to keep the scheme stable, lateral subsurface flow flux/removal,
// and theta(psi,K) with optional soil properties look up table (LUT) & initial hint.
// DEVELOPED IN SUPPORT OF THE DISCRETE SOIL MOISTURE BALANCE MODULE (DSBM) - FLO 9/2025

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "soil_helpers.h"
#include "soil_config.h"

#ifndef THETA_MIN
#define THETA_MIN 1.0e-03   // m3/m3, fictitious minimum water content needed because CH assumes theta_r = 0
#endif


/************/ // Allocate and build a CH-style LUT over Theta = (theta-theta_r)/(theta_sat-theta_r)
int soil_build_ch_lut(
    int n_points,
    double Theta_min_in,               // 0<Theta_min<=1; we clamp internally
    double theta_r,
    double theta_sat,
    double b_exp,
    double phi_sat_cm,
    double K_sat_cm_per_h,
    SoilLookupTables *lut_out)
{
    if (!lut_out) return 1;
    if (n_points < 2) n_points = 2;

    double Theta_min = Theta_min_in;
    if (Theta_min < 1.0e-6) Theta_min = 1.0e-6;
    if (Theta_min > 1.0)    Theta_min = 1.0;

    double phi_sat_m     = phi_sat_cm / 100.0;
    double K_sat_m_per_h = K_sat_cm_per_h / 100.0;

    double lnTheta_min = log(Theta_min);
    double lnTheta_max = 0.0; // ln(1)
    double dln = (lnTheta_max - lnTheta_min) / (double)(n_points - 1);
    if (dln <= 0.0) dln = 1.0; // safety

    double slope_lnpsi = -b_exp;            // ln(psi) = ln(phi_sat) + (-b)*ln(Theta)
    double slope_lnK   =  2.0*b_exp + 3.0;  // ln(K)   = ln(K_sat)   + (2b+3)*ln(Theta)

    double *lnpsi = (double*)malloc((size_t)n_points * sizeof(double));
    double *lnK   = (double*)malloc((size_t)n_points * sizeof(double));
    if (!lnpsi || !lnK) {
        free(lnpsi); free(lnK);
        return 2;
    }

    for (int i = 0; i < n_points; i++) {
        double lnT = lnTheta_min + i * dln;
        lnpsi[i] = log(phi_sat_m)     + slope_lnpsi * lnT;
        lnK[i]   = log(K_sat_m_per_h) + slope_lnK   * lnT;
    }

    lut_out->n = n_points;
    lut_out->lnTheta_min = lnTheta_min;
    lut_out->dlnTheta = dln;
    lut_out->inv_dlnTheta = 1.0 / dln;
    lut_out->lnpsi = lnpsi;
    lut_out->lnK   = lnK;

    (void)theta_r; (void)theta_sat; // not stored; scaling is in Theta
    return 0;
}

/************/ // Free a LUT
void soil_free_ch_lut(SoilLookupTables *lut)
{
    if (!lut) return;
    free(lut->lnpsi);
    free(lut->lnK);
    memset(lut, 0, sizeof(*lut));
}

/************/ // CH analytic relations (psi in m, K in m/h)
double psi_from_theta(double theta,
                      double theta_sat,
                      double phi_sat_cm,
                      double b_exp)
{
    double th = theta;
    if (th < THETA_MIN) th = THETA_MIN;
    if (th > theta_sat) th = theta_sat;

    double phi_m = phi_sat_cm / 100.0;
    double ratio = th / theta_sat;
    double psi_m = phi_m * pow(ratio, -b_exp);
    return psi_m;
}

/************/
double theta_from_psi(double psi_m,
                      double theta_sat,
                      double phi_sat_cm,
                      double b_exp)
{
    if (psi_m <= 0.0) return theta_sat;
    double phi_m = phi_sat_cm / 100.0;
    double ratio = psi_m / phi_m;
    double th = theta_sat * pow(ratio, -1.0 / b_exp);
    if (th < 0.0) th = 0.0;
    if (th > theta_sat) th = theta_sat;
    return th;
}

/************/
double K_from_theta(double theta,
                    double theta_sat,
                    double K_sat_cm_per_h,
                    double b_exp)
{
    double K_sat_m_per_h = K_sat_cm_per_h / 100.0;
    double th = theta;
    if (th < THETA_MIN) th = THETA_MIN;
    if (th > theta_sat) th = theta_sat;

    double ratio = th / theta_sat;
    double expo  = 2.0 * b_exp + 3.0;
    double K = K_sat_m_per_h * pow(ratio, expo);
    if (K < 1.0e-16) K = 1.0e-16;
    return K;
}

/************/ // Evaluate psi,K from theta using LUT if provided; updates hint in place
void eval_theta_to_props(const SoilLookupTables *lut,
                         double theta,
                         double theta_r,
                         double theta_sat,
                         double phi_sat_cm,
                         double K_sat_cm_per_h,
                         double b_exp,
                         int *hint,              // in/out; ignored if lut==NULL
                         double *psi_m_out,
                         double *K_m_per_h_out)
{
    if (lut) {
        // Theta scaled by residual, as in your LUT: Theta=(theta-theta_r)/(theta_sat-theta_r)
        double denom = theta_sat - theta_r;
        if (denom <= 0.0) denom = 1.0;
        double Theta = (theta - theta_r) / denom;
        if (Theta < 0.0) Theta = 0.0;
        if (Theta > 1.0) Theta = 1.0;

        double lnT;
        if (Theta <= 0.0) lnT = lut->lnTheta_min;
        else {
            lnT = log(Theta);
            if (lnT < lut->lnTheta_min) lnT = lut->lnTheta_min;
            if (lnT > 0.0) lnT = 0.0;
        }

        double f = (lnT - lut->lnTheta_min) * lut->inv_dlnTheta;
        int i0;
        if (f <= 0.0) {
            i0 = 0;
            *psi_m_out     = exp(lut->lnpsi[0]);
            *K_m_per_h_out = exp(lut->lnK[0]);
        } else {
            double fn = (double)(lut->n - 1);
            if (f >= fn) {
                i0 = lut->n - 2;
                *psi_m_out     = exp(lut->lnpsi[lut->n - 1]);
                *K_m_per_h_out = exp(lut->lnK[lut->n - 1]);
            } else {
                if (hint && *hint >= 0) {
                    i0 = *hint;
                    while (i0 > 0 && f < (double)i0) i0--;
                    while (i0 < lut->n - 2 && f > (double)(i0 + 1)) i0++;
                } else {
                    i0 = (int)f;
                }
                if (i0 < 0) i0 = 0;
                if (i0 > lut->n - 2) i0 = lut->n - 2;
                double t = f - (double)i0;
                double lnpsi = lut->lnpsi[i0] + t * (lut->lnpsi[i0+1] - lut->lnpsi[i0]);
                double lnK   = lut->lnK[i0]   + t * (lut->lnK[i0+1]   - lut->lnK[i0]);
                *psi_m_out     = exp(lnpsi);
                *K_m_per_h_out = exp(lnK);
            }
        }
        if (hint) *hint = i0;
        (void)phi_sat_cm; (void)K_sat_cm_per_h; (void)b_exp; // not needed in LUT path
    } else {
        *psi_m_out     = psi_from_theta(theta, theta_sat, phi_sat_cm, b_exp);
        *K_m_per_h_out = K_from_theta(theta, theta_sat, K_sat_cm_per_h, b_exp);
        (void)theta_r;
    }
}

/************/ // Vector version for NDISC discs
void compute_props_with_option_stateless(const SoilLookupTables *lut,
                                         const double theta[NDISC],
                                         double psi_m_out[NDISC],
                                         double K_m_per_h_out[NDISC],
                                         double theta_r,
                                         double theta_sat,
                                         double K_sat_cm_per_h,
                                         double phi_sat_cm,
                                         double b_exp,
                                         int hint_inout[NDISC])
{
    for (int i = 0; i < NDISC; i++) {
        int *hptr = hint_inout ? &hint_inout[i] : NULL;
        eval_theta_to_props(lut, theta[i], theta_r, theta_sat,
                            phi_sat_cm, K_sat_cm_per_h, b_exp,
                            hptr, &psi_m_out[i], &K_m_per_h_out[i]);
    }
}

/************/ // Darcy-Buckingham flux; downward positive; m/h
double flux_DB_pair(double psi_up, double K_up,
                    double psi_dn, double K_dn,
                    double dz_up,  double dz_dn)
{
    double dz_int = 0.5 * (dz_up + dz_dn);
    double term   = 1.0 + (psi_dn - psi_up) / dz_int;
    double K_int  = 0.5 * (K_up + K_dn);
    if (K_int < 1.0e-16) K_int = 1.0e-16;
    return K_int * term;
}

/************/ // Choose number of substeps (uses dt_hours)
// Heuristic stability criteria for forward-Euler routing -FLO:
// 1) flux_ratio: how much water would move in dt_hours? 10% of thinnest disc used.
// 2) rain_ratio: fraction of top-cell storage capacity asked for this hour. 20% used.
// The worse of the two determines the substep count (1..12).
int choose_n_sub_dt(double rain_mm_per_h,
                    double theta1, double theta_sat, double dz1,
                    double q12_0, double q23_0, double q34_0,
                    double dzmin,
                    double dt_hours)
{
    int n_sub = 1;

    double qmax0 = fabs(q12_0);
    if (fabs(q23_0) > qmax0) qmax0 = fabs(q23_0);
    if (fabs(q34_0) > qmax0) qmax0 = fabs(q34_0);

    double move_potential = qmax0 * dt_hours;  // m
    double flux_ratio = (dzmin > 0.0) ? (move_potential / (0.10 * dzmin)) : 0.0;

    double rain_rate_m_per_h = rain_mm_per_h / 1000.0;
    double rain_hour_m = rain_rate_m_per_h * dt_hours;
    double cap1_m = (theta_sat - theta1) * dz1;
    if (cap1_m < 1e-12) cap1_m = 1e-12;
    double rain_ratio = rain_hour_m / (0.20 * cap1_m);

    double severity = 0.0;
    if (flux_ratio > severity) severity = flux_ratio;
    if (rain_ratio > severity) severity = rain_ratio;

    if (severity <= 1.0) {
        n_sub = (rain_mm_per_h > 0.0) ? 2 : 1;
    } else if (severity <= 2.0) {
        n_sub = 4;
    } else if (severity <= 3.0) {
        n_sub = 6;
    } else if (severity <= 5.0) {
        n_sub = 8;
    } else {
        n_sub = 12;
    }

    if (n_sub < 1)  n_sub = 1;
    if (n_sub > 12) n_sub = 12;
    return n_sub;
}

/************/ // Lateral removal substep (m removed this substep)
double remove_lateral_to_subsurface_nash_substep(double theta[NDISC], const double dz[NDISC],
                                                 double theta_fc, double theta_sat,
                                                 double rate_const_m_per_h, double dt_sub,
                                                 double removed_by_disc_accum_m[NDISC])
{
    double total_removed_m = 0.0;

    for (int i = 0; i < NDISC; i++) {
        if (theta[i] <= theta_fc) continue;

        double denom = theta_sat - theta_fc;
        if (denom < 1.0e-12) denom = 1.0e-12;

        double frac = (theta[i] - theta_fc) / denom;
        if (frac < 0.0) frac = 0.0;
        if (frac > 1.0) frac = 1.0;

        double potential_m = rate_const_m_per_h * frac * dt_sub;
        double avail_m = (theta[i] - theta_fc) * dz[i];
        if (avail_m < 0.0) avail_m = 0.0;

        double take_m = potential_m;
        if (take_m > avail_m) take_m = avail_m;

        if (take_m > 0.0) {
            theta[i] -= take_m / dz[i];
            if (theta[i] < theta_fc) theta[i] = theta_fc;

            total_removed_m += take_m;
            if (removed_by_disc_accum_m) removed_by_disc_accum_m[i] += take_m;
        }
    }

    return total_removed_m;
}

/************/ // Hydrostatic storage integral (helper)
static double storage_given_zwt(double z_wt,
                                double soil_depth_m,
                                double theta_sat,
                                double phi_sat_cm,
                                double b_exp)
{
    const int    NINT = 4000;
    const double dz   = soil_depth_m / (double)NINT;
    const double phi_sat_m = phi_sat_cm / 100.0;

    double total = 0.0;
    for (int i = 0; i < NINT; i++) {
        double z_mid = (i + 0.5) * dz;
        double z_cf_top = z_wt - phi_sat_m;
        double theta_here;
        if (z_mid >= z_cf_top && z_mid <= z_wt) theta_here = theta_sat;  // capillary fringe saturated
        else if (z_mid > z_wt)                   theta_here = theta_sat;  // below WT
        else {
            double psi_m = z_wt - z_mid; if (psi_m < 0.0) psi_m = 0.0;
            theta_here = theta_from_psi(psi_m, theta_sat, phi_sat_cm, b_exp);
        }
        total += theta_here * dz;
    }
    return total;
}

/************/ // Initialize theta profile from a target total storage [m].
// Uses bisection on water-table depth (zwt) to find the hydrostatic
// equilibrium profile whose integrated storage matches target_storage_m.
// The capillary fringe above the water table is fully saturated; above
// that, theta follows the Clapp-Hornberger retention curve.
void initialize_hydrostatic_from_storage(double soil_depth_m,
                                         double theta_sat,
                                         double phi_sat_cm,
                                         double b_exp,
                                         const double zc[NDISC],
                                         double target_storage_m,
                                         double *zwt_out,
                                         double theta_out[NDISC])
{
    const double phi_sat_m = phi_sat_cm / 100.0;

    double Smin = 0.0;
    double Smax = theta_sat * soil_depth_m;
    if (target_storage_m < Smin) target_storage_m = Smin;
    if (target_storage_m > Smax) target_storage_m = Smax;

    double z_lo = -10.0;
    double z_hi = soil_depth_m + 200.0;

    double S_lo = storage_given_zwt(z_lo, soil_depth_m, theta_sat, phi_sat_cm, b_exp);
    double S_hi = storage_given_zwt(z_hi, soil_depth_m, theta_sat, phi_sat_cm, b_exp);

    if (S_lo < S_hi) {
        double tmpS = S_lo; S_lo = S_hi; S_hi = tmpS;
        double tmpZ = z_lo; z_lo = z_hi; z_hi = tmpZ;
    }

    if (target_storage_m >= S_lo) {
        *zwt_out = z_lo;
    } else if (target_storage_m <= S_hi) {
        *zwt_out = z_hi;
    } else {
        for (int iter = 0; iter < 100; iter++) {
            double z_mid = 0.5 * (z_lo + z_hi);
            double S_mid = storage_given_zwt(z_mid, soil_depth_m, theta_sat, phi_sat_cm, b_exp);
            if (S_mid > target_storage_m) { z_lo = z_mid; S_lo = S_mid; }
            else { z_hi = z_mid; S_hi = S_mid; }
            double width = fabs(z_hi - z_lo);
            double Serr  = fabs(S_mid - target_storage_m);
            if (width < 1.0e-10 || Serr < 1.0e-10) { *zwt_out = z_mid; break; }
            if (iter == 99) *zwt_out = z_mid;
        }
    }

    double z_wt = *zwt_out;
    for (int i = 0; i < NDISC; i++) {
        double z = zc[i];
        double z_cf_top = z_wt - phi_sat_m;
        if (z >= z_cf_top) theta_out[i] = theta_sat;
        else {
            double psi_m = z_wt - z; if (psi_m < 0.0) psi_m = 0.0;
            theta_out[i] = theta_from_psi(psi_m, theta_sat, phi_sat_cm, b_exp);
        }
    }
}
