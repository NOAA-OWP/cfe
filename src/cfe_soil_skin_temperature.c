/*
 * NOAA-OWP/cfe - Version 3 of Conceptual Functional Equivalent to the stormflow/runoff
 *                generation components of the NOAA/NWS National Water Model version 3.1
 *                and earlier
 *
 * Originally conceived and developed by:
 *         Fred L. Ogden, Chief Scientist, NOAA/NWS
 *         Office of Water Prediction, Tuscaloosa, AL
 */

#include <math.h>
#include <stddef.h>

#include "cfe_soil_skin_temperature.h"

#define CFE_PI 3.14159265358979323846

void initialize_soil_temperature_state(
    cfe_pet_temperature_state_struct *temperature_state)
{
    if (temperature_state == NULL) return;

    temperature_state->skin_temperature_k = 279.15;
    temperature_state->upper_soil_temperature_k = 279.15;
    temperature_state->estimated_annual_air_temperature_k = 279.15;
    temperature_state->air_temperature_time_integral_k_s = 0.0;
    temperature_state->accumulated_time_s = 0.0;
    temperature_state->initialized = 1;
}

/*
 * True if dlwrf_surface_w_per_m2 looks like a genuine, physically
 * plausible observed/modeled value rather than the uninitialized
 * sentinel, an unset 0.0 default, or garbage.  See the constants'
 * definitions in cfe_types.h for the reasoning behind the bounds.
 * Fred, 2026-08-06.
 */
static int dlwrf_surface_is_valid(double dlwrf_surface_w_per_m2)
{
    return isfinite(dlwrf_surface_w_per_m2) &&
           dlwrf_surface_w_per_m2 > CFE_DLWRF_SURFACE_MIN_VALID_W_PER_M2 &&
           dlwrf_surface_w_per_m2 < CFE_DLWRF_SURFACE_MAX_VALID_W_PER_M2;
}

static void update_annual_mean_and_upper_soil_temperature(
    double air_temperature_k,
    int day_of_year,
    double timestep_s,
    cfe_pet_temperature_state_struct *temperature_state)
{
    const double initial_annual_mean_k = 279.15;
    const double seconds_per_year = 365.0 * 86400.0;
    double cumulative_mean_air_temperature_k;
    double forcing_weight;
    double seasonal_coupling;

    temperature_state->air_temperature_time_integral_k_s +=
        air_temperature_k * timestep_s;
    temperature_state->accumulated_time_s += timestep_s;

    cumulative_mean_air_temperature_k =
        temperature_state->air_temperature_time_integral_k_s /
        temperature_state->accumulated_time_s;

    forcing_weight = CLAMP(
        temperature_state->accumulated_time_s / seconds_per_year,
        0.0,
        1.0);

    temperature_state->estimated_annual_air_temperature_k =
        (1.0 - forcing_weight) * initial_annual_mean_k +
        forcing_weight * cumulative_mean_air_temperature_k;

    seasonal_coupling =
        0.5 + 0.5 * cos(
            2.0 * CFE_PI * ((double)day_of_year - 172.0) / 365.0);

    temperature_state->upper_soil_temperature_k =
        temperature_state->estimated_annual_air_temperature_k +
        (air_temperature_k -
         temperature_state->estimated_annual_air_temperature_k) *
        seasonal_coupling;
}

static void update_soil_skin_temperature(
    double air_temperature_k,
    double downwelling_shortwave_w_per_m2,
    double dlwrf_surface_w_per_m2,
    double timestep_s,
    cfe_pet_temperature_state_struct *temperature_state)
{
    const double stefan_boltzmann = 5.670374419e-8;
    const double surface_albedo = 0.25;
    const double surface_emissivity = 0.96;
    const double atmospheric_emissivity = 0.757;
    const double air_density_kg_per_m3 = 1.225;
    const double air_specific_heat_j_per_kg_k = 1004.0;
    /* Aerodynamic resistance: CFE_AERODYNAMIC_RESISTANCE_S_PER_M
     * (cfe_types.h), shared with calculate_bare_soil_evap.c. */
    const double soil_thermal_conductivity_w_per_m_k = 1.0;
    const double conduction_distance_m = 0.05;
    const double skin_heat_capacity_j_per_m2_k = 100000.0;

    double old_skin_temperature_k;
    double downwelling_longwave_w_per_m2;
    double radiative_transfer_coefficient_w_per_m2_k;
    double sensible_transfer_coefficient_w_per_m2_k;
    double conductive_transfer_coefficient_w_per_m2_k;
    double storage_coefficient_w_per_m2_k;
    double forcing_term_w_per_m2;
    double denominator_w_per_m2_k;
    double updated_skin_temperature_k;

    old_skin_temperature_k = temperature_state->skin_temperature_k;

    /*
     * Prefer the real forcing when it's available and physically
     * plausible; otherwise fall back to the clear-sky estimate from air
     * temperature.  Fred, 2026-08-06.
     */
    if (dlwrf_surface_is_valid(dlwrf_surface_w_per_m2)) {
        downwelling_longwave_w_per_m2 = dlwrf_surface_w_per_m2;
    }
    else {
        downwelling_longwave_w_per_m2 =
            atmospheric_emissivity * stefan_boltzmann *
            pow(air_temperature_k, 4.0);
    }

    radiative_transfer_coefficient_w_per_m2_k =
        4.0 * surface_emissivity * stefan_boltzmann *
        pow(old_skin_temperature_k, 3.0);

    sensible_transfer_coefficient_w_per_m2_k =
        air_density_kg_per_m3 * air_specific_heat_j_per_kg_k /
        CFE_AERODYNAMIC_RESISTANCE_S_PER_M;

    conductive_transfer_coefficient_w_per_m2_k =
        soil_thermal_conductivity_w_per_m_k /
        conduction_distance_m;

    storage_coefficient_w_per_m2_k =
        skin_heat_capacity_j_per_m2_k / timestep_s;

    forcing_term_w_per_m2 =
        (1.0 - surface_albedo) * downwelling_shortwave_w_per_m2 +
        surface_emissivity * downwelling_longwave_w_per_m2 -
        surface_emissivity * stefan_boltzmann *
            pow(old_skin_temperature_k, 4.0) +
        radiative_transfer_coefficient_w_per_m2_k *
            old_skin_temperature_k +
        sensible_transfer_coefficient_w_per_m2_k *
            air_temperature_k +
        conductive_transfer_coefficient_w_per_m2_k *
            temperature_state->upper_soil_temperature_k;

    denominator_w_per_m2_k =
        storage_coefficient_w_per_m2_k +
        radiative_transfer_coefficient_w_per_m2_k +
        sensible_transfer_coefficient_w_per_m2_k +
        conductive_transfer_coefficient_w_per_m2_k;

    if (denominator_w_per_m2_k <= 0.0 ||
        !isfinite(denominator_w_per_m2_k)) {
        return;
    }

    updated_skin_temperature_k =
        (storage_coefficient_w_per_m2_k * old_skin_temperature_k +
         forcing_term_w_per_m2) /
        denominator_w_per_m2_k;

    temperature_state->skin_temperature_k = CLAMP(
        updated_skin_temperature_k,
        223.15,
        343.15);
}

int update_soil_skin_temperature_state(
    const cfe_forcing_struct *forcing,
    int dt_seconds,
    int day_of_year,
    cfe_state_struct *state)
{
    double air_temperature_k;
    cfe_pet_temperature_state_struct *temperature_state;

    if (forcing == NULL || state == NULL ||
        dt_seconds <= 0 || day_of_year < 1 || day_of_year > 366) {
        return -1;
    }

    air_temperature_k = forcing->TMP_2maboveground;
    if (!isfinite(air_temperature_k) ||
        !isfinite(forcing->DSWRF_surface)) {
        return -1;
    }

    temperature_state = &state->pet_temperature_state;

    if (!temperature_state->initialized) {
        initialize_soil_temperature_state(temperature_state);
        temperature_state->skin_temperature_k = air_temperature_k;
        temperature_state->upper_soil_temperature_k = air_temperature_k;
    }

    update_annual_mean_and_upper_soil_temperature(
        air_temperature_k,
        day_of_year,
        (double)dt_seconds,
        temperature_state);

    update_soil_skin_temperature(
        air_temperature_k,
        forcing->DSWRF_surface,
        forcing->DLWRF_surface,
        (double)dt_seconds,
        temperature_state);

    return 0;
}
