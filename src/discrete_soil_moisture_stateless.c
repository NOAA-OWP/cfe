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
 

// ======================================================================
// Stateless soil step (BMI kernel) NDISC generic
// Positive vertical flux is downward.
// Units:
//   theta: m3/m3
//   psi:   m
//   K:     m/h
//   rainfall, PET inputs: mm/h (converted to m/h inside)
//   step-integrated totals: m
//   rates: m/h
//
// Author: Fred L. Ogden, July, 2025, NOAA/NWS Office of Water Prediction
//
// This code solves the soil moisture evolution given rainfall input with
// percolation output to groundwater plus output to a lateral flow routine
// from a homogeneous 2 m thick soil, discretized into 0.1, 0.3, 0.6, and 
// 1.0 m thick discretizations (discs) (top-down), as in Noah-MP.  It uses
// Darcy-Buckingham flux calculations and the arithmetic average of the
// unsaturated hydraulic conductivity in the discs on either side of the
// interface between them.   Like Noah-MP is uses the field capacity as
// the threshold to activate the percolation and lateral flow fluxes.  It
// extracts AET from the wettest root zone disc.  It uses substeps to keep
// the solution stable, based on a fraction of the available pore space 
// filled by rainfall during a sub-time-step to less than 20%, or the 
// distance that the Darcy flux moves in a sub-time-step to be less than
// 10 percent of the disc thickness. These ratios were determined by trial
// and error, and may not be optimal.  The code can optionally use look-
// up tables of pre-calculated Clapp-Hornberger K(theta) and psi(theta)
// consisting of 5 values (defined in cfe_soil_discrete.h) to reduce 
// computation.  The Clapp-Hornberger functions are linear after logarithmic
// transform so the fit is fantastic, even with only 5 points.
// ======================================================================

#include <math.h>
#include "cfe.h"
#include "discrete_soil_moisture.h"
#include "soil_helpers.h"

#ifndef THETA_MIN
#define THETA_MIN 1.0e-03
#endif

#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);  _a < _b ? _a : _b; })

//####################################
static inline double storage_sum_ndisc(const double *theta, const double *dz_m)
{
    double s = 0.0;
    for (int i = 0; i < NDISC; i++) s += theta[i] * dz_m[i];
    return s;
}

//##############################
static inline int any_disc_above(const double *theta, double thresh, int ndisc)
{
    for (int i = 0; i < ndisc; i++) if (theta[i] > thresh) return 1;
    return 0;
}

// Pick a conservative substep count based on initial fluxes and rainfall demand.
// Generic across NDISC.  Note: discretization is abbreviated here as disc.
//
//
//#############################
static int choose_n_substeps_generic(double dt_hours,
                                double rain_mm_per_h,
                                const double *dz_m,
                                const double *theta, double theta_sat,
                                const double *q0_m_per_h, // [0..NDISC-2]
                                int ndisc)
{
    int nintf = ndisc - 1;  // the number of interfaces between discs

    // Max interface magnitude
    double qmax = 0.0;
    for (int i = 0; i < nintf; i++) {
        double a = fabs(q0_m_per_h[i]);
        if (a > qmax) qmax = a;
    }

    // Thinnest discretization
    double dzmin = dz_m[0];  // Uses Noah-MP discs of 0.1, 0.3, 0.6, 1.0 m

    // Flux criterion: how much of the water would move in dt_hours? 10% used here.
    // These ratios are heuristic stability criteria for forward-Euler routing.  -FLO
    double move_potential = qmax * dt_hours;                 // m
    double flux_ratio     = move_potential / (0.10 * dzmin);

    // Rain criterion: fraction of top-cell storage capacity asked for this hour. 20% used here.
    double rain_rate_m_per_h = rain_mm_per_h / 1000.0;
    double rain_hour_m       = rain_rate_m_per_h * dt_hours;
    double cap1_m            = (theta_sat - theta[0]) * dz_m[0];
    if (cap1_m < 1e-12) cap1_m = 1e-12;
    double rain_ratio        = rain_hour_m / (0.20 * cap1_m);

    double severity = 0.0;
    if (flux_ratio > rain_ratio) { 
        severity = flux_ratio;
    } else {
        severity = rain_ratio;
    }

    int n_substeps = 1;
    if (severity <= 1.0 && rain_mm_per_h > 0.0) n_substeps = 2;
    else if (severity <= 2.0)  n_substeps = 4;
    else if (severity <= 3.0)  n_substeps = 6;
    else if (severity <= 5.0)  n_substeps = 8;
    else                       n_substeps = 12;

    // n_substeps is always 1..12 from the if/else chain above; clamps are defensive only
    if (n_substeps < 1)  n_substeps = 1;
    if (n_substeps > 12) n_substeps = 12;
    return n_substeps;
}

//##############################
int DSBM_step_one_hour_stateless(
    const SoilControl        *control,
    const SoilGeometry       *geom,
    const SoilParameters     *params,
    const SoilLookupTables   *lut,                     // may be NULL
    const SoilStateIn        *state_in,
    struct EVAPOTRANSPIRATION_STRUCTURE *evap_struct,  // Pass the actual ET struct
    const SoilForcing        *forcing,
    SoilStateOut             *state_out,
    SoilFluxes               *flux,
    TimestepSoilVolbal       *volbal,
    FILE                     *debug_fptr)              // may be NULL too
{
    (void)debug_fptr;

    // ---- local working state + zero outputs in one pass ----------------------
    double theta[NDISC];
    for (int i = 0; i < NDISC; i++) {
        theta[i] = state_in->theta_in[i];
        flux->AET_by_disc_m[i] = 0.0;
        flux->lateral_by_disc_m[i] = 0.0;
        flux->interface_vol_m[i] = 0.0;
        flux->interface_rate_m_per_h[i] = 0.0;
        state_out->ch_lut_hint_out[i] = state_in->ch_lut_hint_in[i];
    }
    flux->percolation_to_gw_m = 0.0;
    flux->rain_into_soil_m    = 0.0;
    flux->rain_excess_m       = 0.0;
    flux->n_substeps_used          = 0;

    volbal->in_rain_m       = 0.0;
    volbal->excess_m        = 0.0;
    volbal->perc_m          = 0.0;
    volbal->AET_m           = 0.0;
    volbal->lateral_m       = 0.0;
    volbal->delta_storage_m = 0.0;
    volbal->residual_m      = 0.0;

    const int ndisc  = control->ndisc;
    const int nintf  = ndisc - 1;
    const double delta_t_h = control->dt_hours;

    const double rain_rate_m_per_h = forcing->rain_mm_per_h / 1000.0;

    const double theta_floor = fmax(THETA_MIN, params->theta_r);

    // initial storage
    const double storage_start = storage_sum_ndisc(theta, geom->dz_m);

    //-- NEW
    // Create a working copy of the input soil state, and populate theta in discs

    SoilStateIn temp_soil_state = *state_in;

    for (int i = 0; i < NDISC; i++) {
        temp_soil_state.theta_in[i] = theta[i];
    }

    /*
     * Remove area-weighted bare-soil evaporation from disc 1, then apply
     * forest/root-zone AET using the existing wettest-root-disc method.
     * The externally calculated bare-soil flux is capped again here so the
     * state update cannot cross the wilting-point water content.
     */
    double requested_bare_aet_m =
        evap_struct->actual_bare_soil_evaporation_m_per_timestep;
    double available_top_water_m =
        fmax(temp_soil_state.theta_in[0] - params->theta_wp, 0.0) *
        geom->dz_m[0];
    double actual_bare_aet_m =
        fmin(fmax(requested_bare_aet_m, 0.0), available_top_water_m);

    temp_soil_state.theta_in[0] -= actual_bare_aet_m / geom->dz_m[0];
    evap_struct->actual_bare_soil_evaporation_m_per_timestep =
        actual_bare_aet_m;

    et_from_soil_discrete(control, geom, params, &temp_soil_state, evap_struct);

    // Update theta array with post-ET values
    for (int i = 0; i < NDISC; i++) {
        theta[i] = temp_soil_state.theta_in[i];
    }

    // Track ET removal by disc for flux accounting
    double et_removed = 0.0;
    for (int i = 0; i < NDISC; i++) {
        double et_removed_from_disc_m =
            (state_in->theta_in[i] - theta[i]) * geom->dz_m[i];
        et_removed += et_removed_from_disc_m;
        flux->AET_by_disc_m[i] = et_removed_from_disc_m;
    }

    volbal->AET_m = et_removed;
    //-- END NEW

    // Calculate initial fluxes at disc interfaces
    double q0[NDISC-1];  
    for (int i = 0; i < nintf; i++) {
        // calculate the Darcy-Buckingham (DB) flux from disc 0-1, disc 1-2, disc 2-3.
        q0[i] = flux_DB_pair(state_in->psi_in[i], state_in->K_in[i],
                             state_in->psi_in[i+1], state_in->K_in[i+1],
                             geom->dz_m[i], geom->dz_m[i+1]);
    }

    // Determine the number of substeps (needed in case of very wet soils in disc1)
    int n_substeps = choose_n_substeps_generic(delta_t_h, forcing->rain_mm_per_h,
                                     geom->dz_m, theta, params->theta_sat, q0, ndisc);
    // choose function guarantees 1..12, but clamp defensively
    if (n_substeps < 1) n_substeps = 1;
    flux->n_substeps_used = n_substeps;

    const double dt_sub = delta_t_h / (double)n_substeps;

    // per-substep work arrays
    double psi[NDISC], K[NDISC];
    double store_cap[NDISC];
    double pot_downflux[NDISC > 1 ? NDISC-1 : 1];  // >=0
    double V_if[NDISC > 1 ? NDISC - 1 : 1];          // signed desired
    double Accept[NDISC + 1] = {0};                 // [0..ndisc], nd = bottom; C99 {0} zeros all elements

    // External inflow (incident) is available to the caller via forcing and delta_t_h.
    // Here we only track infiltrated vs excess ffor the step.
    // Loop over substeps

    //  <----------------------------------------------------------- Start of substep loop
    for (int substep = 0; substep < n_substeps; substep++) {
        // See if partitioned soil moisture from Schaake/Xinanjiang fits into disc 0
        double rain_sub = rain_rate_m_per_h * dt_sub;  // m
        double cap1     = (params->theta_sat - theta[0]) * geom->dz_m[0];
        if (cap1 < 0.0) cap1 = 0.0;

        double used   = rain_sub;
        if (used > cap1) used = cap1;

        double excess = rain_sub - used;
        if (excess < 0.0) excess = 0.0;

        theta[0] += used / geom->dz_m[0];
        if (theta[0] > params->theta_sat) theta[0] = params->theta_sat;

        flux->rain_into_soil_m += used;
        flux->rain_excess_m    += excess;

        // Properties after rainfall addition (LUT or analytic)
        int hint_local[NDISC];
        for (int i = 0; i < NDISC; i++) hint_local[i] = state_out->ch_lut_hint_out[i];

        compute_props_with_option_stateless(
            (control->use_ch_lookup_table ? lut : NULL),
            theta, psi, K,
            params->theta_r, params->theta_sat,
            params->K_sat_cm_per_h, params->phi_sat_cm, params->b_exp,
            hint_local);

        // keep updated hints
        for (int i = 0; i < NDISC; i++) state_out->ch_lut_hint_out[i] = hint_local[i];

        // Interface fluxes and desired substep volumes
        for (int i = 0; i < nintf; i++) {
            const double q = flux_DB_pair(psi[i], K[i], psi[i+1], K[i+1],
                                          geom->dz_m[i], geom->dz_m[i+1]);
            V_if[i] = q * dt_sub;
            pot_downflux[i] = (q > 0.0) ? (q * dt_sub) : 0.0;
        }

        // Per-disc free storage up to saturation
        for (int d = 0; d < ndisc; d++) {
            double s = (params->theta_sat - theta[d]) * geom->dz_m[d];
            if (s < 0.0) s = 0.0;
            store_cap[d] = s;
        }

        // bottom potential percolation (K(theta_bottom) * limiter)
        // Not limited by GW reservoir storage: the exponential/nonlinear
        // reservoir has no true capacity ceiling — storage_max_m is a
        // curve-shape parameter, not a hard bucket size.
        double bottom_potential = 0.0;
        if (theta[ndisc-1] > params->theta_fc) {
            double K_now = K_from_theta(theta[ndisc-1], params->theta_sat,
                                        params->K_sat_cm_per_h, params->b_exp);
            double percolation_rate = params->perc_limiter_0_to_1 * K_now;     // m/h
            if (percolation_rate > 0.0) {
                bottom_potential = percolation_rate * dt_sub;
            }
        }

        // Downstream acceptance (bottom up)
        // Accept[ndisc] = last store + bottom_potential (now GW-storage-limited)
        //
        // Accept[k] tracks how much water disc k can receive and pass downward.
        // Sweeping bottom-up ensures each disc knows its downstream capacity
        // before accepting water from above, preventing overshoot in a substep.
        Accept[ndisc] = store_cap[ndisc-1] + bottom_potential;

        for (int i = nintf-1; i >= 0; i--) {
            double pass = pot_downflux[i];
            int down_index = i + 2;                 // downstream Accept slot; last interface -> ndisc
            if (down_index >= ndisc) down_index = ndisc;
            if (pass > Accept[down_index]) pass = Accept[down_index];
            Accept[i+1] = store_cap[i] + pass;
        }

        // Apply capped transfers across all interfaces
        // Each transfer is limited by donor availability, receiver capacity,
        // and downstream chain acceptance computed above.
        for (int i = 0; i < nintf; i++) {
            double V = V_if[i];

            double accept_down =
                (i + 2 <= ndisc-1) ? Accept[i+2] : Accept[ndisc];

            if (V > 0.0) {
                double donor_avail = (theta[i] - theta_floor) * geom->dz_m[i];
                if (donor_avail < 0.0) donor_avail = 0.0;

                double recv_space = store_cap[i+1];

                double chain_pass = pot_downflux[i];
                if (chain_pass > accept_down) chain_pass = accept_down;

                double max_out = store_cap[i] + chain_pass;

                if (V > max_out) V = max_out;
                if (V > donor_avail) V = donor_avail;
                if (V > recv_space) V = recv_space;
                if (V < 0.0) V = 0.0;

                theta[i]   -= V / geom->dz_m[i];
                theta[i+1] += V / geom->dz_m[i+1];

                if (theta[i]   < theta_floor)      theta[i]   = theta_floor;
                if (theta[i+1] > params->theta_sat)   theta[i+1] = params->theta_sat;
            }
            else if (V < 0.0) {
                double need = -V;

                double donor_avail = (theta[i+1] - theta_floor) * geom->dz_m[i+1];
                if (donor_avail < 0.0) donor_avail = 0.0;

                double recv_space = store_cap[i];

                double move = need;
                if (move > donor_avail) move = donor_avail;
                if (move > recv_space)  move = recv_space;

                V = -move;

                theta[i+1] -= move / geom->dz_m[i+1];
                theta[i]   += move / geom->dz_m[i];

                if (theta[i+1] < theta_floor)      theta[i+1] = theta_floor;
                if (theta[i]   > params->theta_sat)   theta[i]   = params->theta_sat;
            }

            // accumulate internal interface volume
            flux->interface_vol_m[i] += V;
        }

        // Apply bottom percolation
        double perc_vol = 0.0;
        if (theta[ndisc-1] > params->theta_fc && bottom_potential > 0.0) {
            double avail = (theta[ndisc-1] - params->theta_fc) * geom->dz_m[ndisc-1];
            if (avail < 0.0) avail = 0.0;

            perc_vol = bottom_potential;
            if (perc_vol > avail) perc_vol = avail;

            if (perc_vol > 0.0) {
                theta[ndisc-1] -= perc_vol / geom->dz_m[ndisc-1];
                if (theta[ndisc-1] < params->theta_fc) theta[ndisc-1] = params->theta_fc;
            }
        }
        flux->percolation_to_gw_m += perc_vol;
        flux->interface_vol_m[ndisc-1] += perc_vol;


        // Remove calculated lateral flow to subsurface Nash cascade but only iff any disc > FC
        if (params->klf_per_h > 0.0 && any_disc_above(theta, params->theta_fc, ndisc)) {
            double lat_removed =
                remove_lateral_to_subsurface_nash_substep(
                    theta, geom->dz_m, params->theta_fc, params->theta_sat,
                    params->klf_per_h, dt_sub, flux->lateral_by_disc_m);

            volbal->lateral_m += lat_removed;
        }


        // Safety clamp: theta must stay in [theta_floor, theta_sat].
        // NOTE: clamping introduces a small mass imbalance (water created or
        // destroyed) that appears in the timestep residual. In practice this is
        // negligible because the acceptance logic above prevents large overshoots.
        for (int i = 0; i < ndisc; i++) {
            if (theta[i] < theta_floor)       theta[i] = theta_floor;
            if (theta[i] > params->theta_sat) theta[i] = params->theta_sat;
        }
    } // <------------------------------------------------------------------end substep loop

    // Finalize: write state out and compute interface rates in one pass
    for (int i = 0; i < NDISC; i++) {
        state_out->theta_out[i] = theta[i];
        flux->interface_rate_m_per_h[i] = flux->interface_vol_m[i] / delta_t_h;
    }
    state_out->total_storage_m = storage_sum_ndisc(theta, geom->dz_m);
    state_out->storage_deficit_m = geom->depth_m * params->theta_sat - state_out->total_storage_m;
    
    // volume balances
    double storage_end = storage_sum_ndisc(theta, geom->dz_m);

    volbal->in_rain_m  = flux->rain_into_soil_m;
    volbal->excess_m   = flux->rain_excess_m;
    volbal->perc_m     = flux->percolation_to_gw_m;

    // total lateral already accumulated in volbal->lateral_m inside loop
    // Ensure it matches sum of per-disc laterals (defensive)
    double lat_sum = 0.0; 
    for (int i = 0; i < NDISC; i++) lat_sum += flux->lateral_by_disc_m[i];
    volbal->lateral_m = lat_sum;

    volbal->delta_storage_m = (storage_end - storage_start);

    volbal->residual_m = (volbal->in_rain_m)
                   - (volbal->perc_m + volbal->AET_m + volbal->lateral_m)
                   -  volbal->delta_storage_m;

    return 0;
}


//##############################################################
//############   ET FROM SOIL DISCRETE  ########################
//##############################################################
void et_from_soil_discrete
    (
    const SoilControl*    soil_control,
    const SoilGeometry*   soil_geometry,
    const SoilParameters* soil_parameters,
    SoilStateIn*    soil_state,
    struct EVAPOTRANSPIRATION_STRUCTURE* evap_struct
    )
{
    // FLO September, 2026
    // Distribute PET equally among the root-zone discretizations (discs),
    // consistent with the DSBM formulation used to mimic Noah-MP soil moisture.
    // Each root-zone disc has its own soil-moisture stress multiplier:
    //     transpiration/PET = 0 when theta <= wilting point,
    //     transpiration/PET = 1 when theta >= field capacity,
    //     transpiration/PET varies linearly between wilting point and field capacity.
    // Extraction from any disc is limited so theta cannot fall below the
    // wilting-point moisture content. Unmet demand from one disc is not
    // reassigned to another disc.

    int nroot = soil_control->deepest_root_disc;
    if (nroot < 1) nroot = 1;
    if (nroot > NDISC) nroot = NDISC;

    const double PET = evap_struct->reduced_potential_et_m_per_timestep;
    const double root_fraction = 1.0 / (double)nroot;
    const double stress_denominator =
        fmax(soil_parameters->theta_fc - soil_parameters->theta_wp, 1.0e-12);

    double actual_transpiration_m = 0.0;

    for (int i = 0; i < nroot; i++) {
        const double theta = soil_state->theta_in[i];
        double moisture_stress_multiplier;

        if (theta <= soil_parameters->theta_wp) {
            moisture_stress_multiplier = 0.0;
        } else if (theta >= soil_parameters->theta_fc) {
            moisture_stress_multiplier = 1.0;
        } else {
            moisture_stress_multiplier =
                (theta - soil_parameters->theta_wp) / stress_denominator;
        }

        const double demand_m =
            PET * root_fraction * moisture_stress_multiplier;
        const double available_m =
            fmax(theta - soil_parameters->theta_wp, 0.0) * soil_geometry->dz_m[i];
        const double actual_transpiration_from_disc_m = fmin(demand_m, available_m);

        if (actual_transpiration_from_disc_m > 0.0) {
            soil_state->theta_in[i] -= actual_transpiration_from_disc_m / soil_geometry->dz_m[i];
            if (soil_state->theta_in[i] < soil_parameters->theta_wp)
                soil_state->theta_in[i] = soil_parameters->theta_wp;
        }

        actual_transpiration_m += actual_transpiration_from_disc_m;
    }

    evap_struct->actual_et_from_soil_m_per_timestep = actual_transpiration_m;
    return;
}

