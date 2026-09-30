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

#include "cfe_pet_priestley_taylor.h"

/*
 * True if dlwrf_surface_w_per_m2 looks like a genuine, physically
 * plausible observed/modeled value rather than the uninitialized
 * sentinel, an unset 0.0 default, or garbage.  Kept in sync with the
 * identical check in cfe_soil_skin_temperature.c; see the constants'
 * definitions in cfe_types.h for the reasoning behind the bounds.
 * Fred, 2026-08-06.
 */
static int dlwrf_surface_is_valid(double dlwrf_surface_w_per_m2)
{
    return isfinite(dlwrf_surface_w_per_m2) &&
           dlwrf_surface_w_per_m2 > CFE_DLWRF_SURFACE_MIN_VALID_W_PER_M2 &&
           dlwrf_surface_w_per_m2 < CFE_DLWRF_SURFACE_MAX_VALID_W_PER_M2;
}

double calculate_pet_priestley_taylor(
    const cfe_forcing_struct *forcing,
    int dt_seconds,
    double alpha_pt,
    const cfe_state_struct *state)
{
    const double stefan_boltzmann = 5.670374419e-8;
    const double surface_albedo = 0.25;
    const double surface_emissivity = 0.96;
    const double atmospheric_emissivity = 0.757;
    const double latent_heat_vaporization_j_per_kg = 2.45e6;
    const double psychrometric_constant_kpa_per_k = 0.0665;
    const double liquid_water_density_kg_per_m3 = 998.0;
    const double soil_thermal_conductivity_w_per_m_k = 1.0;
    const double conduction_distance_m = 0.05;
    double air_temperature_k;
    double air_temperature_c;
    double saturation_vapor_pressure_kpa;
    double saturation_vapor_pressure_slope_kpa_per_k;
    double downwelling_longwave_w_per_m2;
    double outgoing_longwave_w_per_m2;
    double net_radiation_w_per_m2;
    double ground_heat_flux_w_per_m2;
    double available_energy_w_per_m2;
    double pet_rate_m_per_s;
    double pet_m_per_timestep;
    double maximum_pet_m_per_timestep;

    const cfe_pet_temperature_state_struct *temperature_state;

    if (forcing == NULL || state == NULL ||
        dt_seconds <= 0 ||
        !isfinite(alpha_pt) || alpha_pt <= 0.0) {
        return 0.0;
    }

    temperature_state = &state->pet_temperature_state;

    air_temperature_k = forcing->TMP_2maboveground;
    if (!isfinite(air_temperature_k) ||
        !isfinite(forcing->DSWRF_surface)) {
        return 0.0;
    }

    if (!temperature_state->initialized) return 0.0;

    downwelling_longwave_w_per_m2 =
        dlwrf_surface_is_valid(forcing->DLWRF_surface)
            ? forcing->DLWRF_surface
            : atmospheric_emissivity * stefan_boltzmann *
                  pow(air_temperature_k, 4.0);

    outgoing_longwave_w_per_m2 =
        surface_emissivity * stefan_boltzmann *
        pow(temperature_state->skin_temperature_k, 4.0);

    net_radiation_w_per_m2 =
        (1.0 - surface_albedo) * forcing->DSWRF_surface +
        surface_emissivity * downwelling_longwave_w_per_m2 -
        outgoing_longwave_w_per_m2;

    ground_heat_flux_w_per_m2 =
        soil_thermal_conductivity_w_per_m_k *
        (temperature_state->skin_temperature_k -
         temperature_state->upper_soil_temperature_k) /
        conduction_distance_m;

    available_energy_w_per_m2 =
        net_radiation_w_per_m2 - ground_heat_flux_w_per_m2;

    if (available_energy_w_per_m2 <= 0.0) return 0.0;

    air_temperature_c = air_temperature_k - 273.15;
    saturation_vapor_pressure_kpa =
        0.6108 * exp(
            (17.27 * air_temperature_c) /
            (air_temperature_c + 237.3));

    saturation_vapor_pressure_slope_kpa_per_k =
        4098.0 * saturation_vapor_pressure_kpa /
        pow(air_temperature_c + 237.3, 2.0);

    pet_rate_m_per_s =
        alpha_pt *
        (saturation_vapor_pressure_slope_kpa_per_k /
         (saturation_vapor_pressure_slope_kpa_per_k +
          psychrometric_constant_kpa_per_k)) *
        (available_energy_w_per_m2 /
         latent_heat_vaporization_j_per_kg) /
        liquid_water_density_kg_per_m3;

    pet_m_per_timestep = pet_rate_m_per_s * (double)dt_seconds;

    maximum_pet_m_per_timestep =
        0.01 * ((double)dt_seconds / 43200.0);

    if (pet_m_per_timestep > maximum_pet_m_per_timestep) {
        pet_m_per_timestep = maximum_pet_m_per_timestep;
    }

    if (pet_m_per_timestep < 0.0 || !isfinite(pet_m_per_timestep)) {
        return 0.0;
    }

    return pet_m_per_timestep;
}
