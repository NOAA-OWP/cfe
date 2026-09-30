/*
 * NOAA-OWP/cfe - Version 3 of Conceptual Functional Equivalent to the stormflow/runoff
 *                generation components of the NOAA/NWS National Water Model version 3.1
 *                and earlier
 *
 * Originally conceived and developed by:
 *         Fred L. Ogden, Chief Scientist, NOAA/NWS
 *         Office of Water Prediction, Tuscaloosa, AL
 *
 * Bare-soil actual evaporation using the Sakaguchi-Zeng soil-surface
 * resistance formulation used as an option in Noah-MP.
 *
 * The function returns evaporation as a positive water depth in meters per
 * model timestep.  It does not modify soil moisture; the caller must remove
 * the returned depth from the upper soil disc and include it in the CFE water
 * balance.
 *
 * Method summary
 * --------------
 * 1. Use liquid volumetric water content in the upper soil disc to estimate
 *    the thickness of the dry surface layer.
 * 2. Divide that dry-layer thickness by an effective vapor diffusivity to
 *    obtain the Sakaguchi-Zeng surface resistance, r_surf [s/m].
 * 3. Reduce saturation vapor pressure at the soil surface with the Kelvin
 *    relation, using Clapp-Hornberger matric potential to obtain RHSUR.
 * 4. Compute vapor transfer through the series resistance r_a + r_surf:
 *
 *        E = rho_air * (q_surface - q_air) / (r_a + r_surf)
 *
 * 5. Convert E from kg/m2/s to m/timestep and limit extraction to liquid
 *    water stored above the wilting-point content in the upper soil disc.
 *
 * Important assumptions
 * ---------------------
 * - The upper CFE soil disc represents the evaporating soil reservoir.
 * - Soil ice is not handled explicitly here; the caller should provide liquid
 *   water content or suppress evaporation when the surface is frozen.
 * - Aerodynamic resistance is supplied by the caller.  This keeps atmospheric
 *   transfer physics separate from the Sakaguchi-Zeng soil resistance.
 * - RSURF_EXP is supplied by the caller.  Noah-MP commonly treats it as a
 *   vegetation/soil parameter rather than a universal constant.
 * - The wilting point is used as the minimum extractable water content.  This
 *   is a CFE water-balance safeguard, not part of the resistance equation.
 *
 * Development notes
 * -----------------
 * This routine was included in CFE to make it more functionally equivalent to
 * WRF-Hydro as employed in the NWM.  Developed by F.L. Ogden with considerable
 * coding help from OpenAI ChatGPT, July/August, 2026.
 *
 */

#include <math.h>
#include <stddef.h>

#include "calculate_bare_soil_evap.h"

static double saturation_vapor_pressure_pa(double temperature_k)
{
    double temperature_c;

    temperature_c = temperature_k - 273.15;

    /* Tetens equation; valid for ordinary land-surface temperatures. */
    return 610.8 * exp((17.27 * temperature_c) /
                       (temperature_c + 237.3));
}

static double specific_humidity_from_vapor_pressure(
    double vapor_pressure_pa,
    double air_pressure_pa)
{
    const double epsilon = 0.622;
    double denominator_pa;

    denominator_pa = air_pressure_pa -
                     (1.0 - epsilon) * vapor_pressure_pa;

    if (denominator_pa <= 0.0) return 0.0;

    return epsilon * vapor_pressure_pa / denominator_pa;
}

/*
 * Calculate bare-soil actual evaporation using Sakaguchi-Zeng resistance.
 *
 * Inputs
 * ------
 * forcing
 *     Meteorological forcing.  Required fields are PRES_surface,
 *     SPFH_2maboveground, and TMP_2maboveground.
 * parameters
 *     CFE soil properties.  Required fields are effective_porosity,
 *     wilting_point, soil_b, and sat_capillary_head_m.
 * state
 *     CFE state.  The upper-disc water content is taken from
 *     soil_discrete_storage_theta[0], the soil skin temperature from
 *     soil skin temperature state, and upper-disc thickness from
 *     soil_geometry.dz_m[0].
 * dt_seconds
 *     Model timestep [s].
 * aerodynamic_resistance_s_per_m
 *     Atmospheric resistance to water-vapor transport [s/m].
 * rsurf_exp
 *     Sakaguchi-Zeng dry-layer shape exponent [-].
 * surface_resistance_s_per_m_out
 *     Optional diagnostic output for r_surf [s/m]; may be NULL.
 * surface_relative_humidity_out
 *     Optional diagnostic output for RHSUR [0-1]; may be NULL.
 *
 * Return value
 * ------------
 * Positive bare-soil evaporation depth [m/timestep].  Returns zero for
 * invalid inputs, condensation conditions, frozen/empty available storage,
 * or nonphysical intermediate values.
 *
 */
double calculate_bare_soil_evap(
    const cfe_forcing_struct *forcing,
    const cfe_parameters_struct *parameters,
    const cfe_state_struct *state,
    int dt_seconds,
    double aerodynamic_resistance_s_per_m,
    double rsurf_exp,
    double *surface_resistance_s_per_m_out,
    double *surface_relative_humidity_out)
{
    const double molecular_diffusivity_water_vapor_m2_per_s = 2.2e-5;
    const double gravitational_acceleration_m_per_s2 = 9.80665;
    const double gas_constant_water_vapor_j_per_kg_k = 461.5;
    const double dry_layer_denominator = 1.71828182845904523536; /* e - 1 */
    const double minimum_effective_saturation = 1.0e-6;
    const double maximum_surface_resistance_s_per_m = 1.0e12;
    const double liquid_water_density_kg_per_m3 = 1000.0;
    const double dry_air_gas_constant_j_per_kg_k = 287.05;

    double air_temperature_k;
    double skin_temperature_k;
    double air_pressure_pa;
    double air_specific_humidity;
    double air_density_kg_per_m3;
    double theta_upper;
    double theta_saturated;
    double theta_wilting;
    double effective_saturation;
    double upper_disc_thickness_m;
    double dry_layer_thickness_m;
    double effective_vapor_diffusivity_m2_per_s;
    double surface_resistance_s_per_m;
    double matric_potential_m;
    double surface_relative_humidity;
    double surface_saturation_vapor_pressure_pa;
    double surface_vapor_pressure_pa;
    double surface_specific_humidity;
    double total_resistance_s_per_m;
    double evaporation_mass_flux_kg_per_m2_s;
    double evaporation_depth_m;
    double available_water_depth_m;
    double diffusivity_exponent;

    if (surface_resistance_s_per_m_out != NULL) {
        *surface_resistance_s_per_m_out = 0.0;
    }
    if (surface_relative_humidity_out != NULL) {
        *surface_relative_humidity_out = 0.0;
    }

    if (forcing == NULL || parameters == NULL || state == NULL ||
        dt_seconds <= 0 ||
        !isfinite(aerodynamic_resistance_s_per_m) ||
        aerodynamic_resistance_s_per_m <= 0.0 ||
        !isfinite(rsurf_exp) || rsurf_exp <= 0.0) {
        return 0.0;
    }

    air_temperature_k = forcing->TMP_2maboveground;
    skin_temperature_k =
        state->pet_temperature_state.skin_temperature_k;
    air_pressure_pa = forcing->PRES_surface;
    air_specific_humidity = forcing->SPFH_2maboveground;

    theta_upper = state->soil_discrete_storage_theta[0];
    theta_saturated = parameters->effective_porosity;
    theta_wilting = parameters->wilting_point;
    upper_disc_thickness_m = state->soil_geometry.dz_m[0];

    if (!isfinite(air_temperature_k) || air_temperature_k <= 0.0 ||
        !isfinite(skin_temperature_k) || skin_temperature_k <= 0.0 ||
        !isfinite(air_pressure_pa) || air_pressure_pa <= 0.0 ||
        !isfinite(air_specific_humidity) || air_specific_humidity < 0.0 ||
        !isfinite(theta_upper) ||
        !isfinite(theta_saturated) || theta_saturated <= 0.0 ||
        !isfinite(theta_wilting) || theta_wilting < 0.0 ||
        theta_wilting >= theta_saturated ||
        !isfinite(parameters->soil_b) || parameters->soil_b <= 0.0 ||
        !isfinite(parameters->sat_capillary_head_m) ||
        parameters->sat_capillary_head_m <= 0.0) {
        return 0.0;
    }

    if (!isfinite(upper_disc_thickness_m) ||
        upper_disc_thickness_m <= 0.0) {
        upper_disc_thickness_m = 0.10;
    }

    effective_saturation = CLAMP(
        theta_upper / theta_saturated,
        minimum_effective_saturation,
        1.0);

    /* Sakaguchi-Zeng dry surface-layer thickness. */
    dry_layer_thickness_m =
        upper_disc_thickness_m *
        (exp(pow(1.0 - effective_saturation, rsurf_exp)) - 1.0) /
        dry_layer_denominator;

    /*
     * Effective vapor diffusivity through the dry soil pore space.  This is
     * the Noah-MP Sakaguchi-Zeng expression written using CFE soil parameter names.
     */
    diffusivity_exponent = 2.0 + 3.0 / parameters->soil_b;
    effective_vapor_diffusivity_m2_per_s =
        molecular_diffusivity_water_vapor_m2_per_s *
        theta_saturated * theta_saturated *
        pow(1.0 - theta_wilting / theta_saturated,
            diffusivity_exponent);

    if (!isfinite(effective_vapor_diffusivity_m2_per_s) ||
        effective_vapor_diffusivity_m2_per_s <= 0.0) {
        return 0.0;
    }

    surface_resistance_s_per_m =
        dry_layer_thickness_m /
        effective_vapor_diffusivity_m2_per_s;

    surface_resistance_s_per_m = CLAMP(
        surface_resistance_s_per_m,
        0.0,
        maximum_surface_resistance_s_per_m);

    /*
     * Clapp-Hornberger matric potential and Kelvin relative humidity.
     * CFE stores saturated capillary head as a positive magnitude, so the
     * matric potential is explicitly negative here.
     */
    matric_potential_m =
        -parameters->sat_capillary_head_m *
        pow(effective_saturation, -parameters->soil_b);

    surface_relative_humidity = exp(
        matric_potential_m * gravitational_acceleration_m_per_s2 /
        (gas_constant_water_vapor_j_per_kg_k * skin_temperature_k));

    surface_relative_humidity = CLAMP(
        surface_relative_humidity,
        0.0,
        1.0);

    surface_saturation_vapor_pressure_pa =
        saturation_vapor_pressure_pa(skin_temperature_k);
    surface_vapor_pressure_pa =
        surface_relative_humidity * surface_saturation_vapor_pressure_pa;

    surface_specific_humidity = specific_humidity_from_vapor_pressure(
        surface_vapor_pressure_pa,
        air_pressure_pa);

    /* Ideal-gas moist-air approximation, adequate for the bulk flux. */
    air_density_kg_per_m3 =
        air_pressure_pa /
        (dry_air_gas_constant_j_per_kg_k * air_temperature_k);

    total_resistance_s_per_m =
        aerodynamic_resistance_s_per_m + surface_resistance_s_per_m;

    if (!isfinite(total_resistance_s_per_m) ||
        total_resistance_s_per_m <= 0.0 ||
        !isfinite(air_density_kg_per_m3) ||
        air_density_kg_per_m3 <= 0.0) {
        return 0.0;
    }

    evaporation_mass_flux_kg_per_m2_s =
        air_density_kg_per_m3 *
        (surface_specific_humidity - air_specific_humidity) /
        total_resistance_s_per_m;

    /* Condensation is not treated as negative AET in this routine. */
    if (!isfinite(evaporation_mass_flux_kg_per_m2_s) ||
        evaporation_mass_flux_kg_per_m2_s <= 0.0) {
        if (surface_resistance_s_per_m_out != NULL) {
            *surface_resistance_s_per_m_out =
                surface_resistance_s_per_m;
        }
        if (surface_relative_humidity_out != NULL) {
            *surface_relative_humidity_out =
                surface_relative_humidity;
        }
        return 0.0;
    }

    evaporation_depth_m =
        evaporation_mass_flux_kg_per_m2_s * (double)dt_seconds /
        liquid_water_density_kg_per_m3;

    available_water_depth_m =
        fmax(theta_upper - theta_wilting, 0.0) *
        upper_disc_thickness_m;

    if (evaporation_depth_m > available_water_depth_m) {
        evaporation_depth_m = available_water_depth_m;
    }

    if (surface_resistance_s_per_m_out != NULL) {
        *surface_resistance_s_per_m_out = surface_resistance_s_per_m;
    }
    if (surface_relative_humidity_out != NULL) {
        *surface_relative_humidity_out = surface_relative_humidity;
    }

    if (!isfinite(evaporation_depth_m) || evaporation_depth_m < 0.0) {
        return 0.0;
    }

    return evaporation_depth_m;
}
