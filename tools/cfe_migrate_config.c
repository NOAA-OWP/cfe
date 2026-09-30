/*
 * cfe_migrate_config — Convert a CFE v2 legacy config to v3 format.
 *
 * Usage: cfe_migrate_config <input_v2_config> <output_v3_config>
 *
 * This is a standalone utility with no dependency on the CFE library.
 * It reads the v2 key=value[units] format, applies key renaming and
 * unit conversions, and writes a v3 config file.
 *
 * Incompatible configs (e.g., Nash surface routing without GIUH) produce
 * an error message and nonzero exit code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAX_LINE 2048
#define MAX_KEY  256
#define MAX_VAL  1024
#define MAX_GIUH 11
#define MAX_NASH 2
#define FMT "%.15e"

/*
 * Estimate wilting point moisture content from Clapp-Hornberger soil parameters.
 * Wilting point is defined as the moisture content at 15 atmospheres of suction.
 *
 *   theta_wp = porosity * (psi_sat / psi_wp)^(1/b)
 *
 * where psi_wp = 15 atm converted to cm of water head.
 * Returns 0.0 if inputs are invalid (caller should use a fallback).
 */
static double wilting_point_from_clapp_hornberger(double porosity, double satpsi_cm, double b_exp)
{
    if (satpsi_cm <= 0.0 || b_exp <= 0.0 || porosity <= 0.0)
        return 0.0;
    double psi_wp_cm = 15.0 * 101325.0 / (9.81 * 998.0) * 100.0;  /* 15 atm in cm H2O */
    return porosity * pow(satpsi_cm / psi_wp_cm, 1.0 / b_exp);
}

/* ---- helpers ---- */

static void trim(char *s)
{
    char *start = s;
    while (isspace((unsigned char)*start)) start++;
    if (start != s) memmove(s, start, strlen(start) + 1);
    char *end = s + strlen(s) - 1;
    while (end >= s && isspace((unsigned char)*end)) *end-- = '\0';
}

static int ci_eq(const char *a, const char *b)
{
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return 0;
        a++; b++;
    }
    return *a == '\0' && *b == '\0';
}

/* Parse "key = value [units]".  Strips comments starting with # or // after value. */
static int parse_line(const char *line, char *key, char *val)
{
    /* skip comments and blank lines */
    const char *p = line;
    while (isspace((unsigned char)*p)) p++;
    if (*p == '#' || *p == '\0' || (*p == '/' && *(p+1) == '/')) return 0;

    char *eq = strchr(p, '=');
    if (!eq) return 0;

    /* key */
    int klen = (int)(eq - p);
    if (klen >= MAX_KEY) klen = MAX_KEY - 1;
    strncpy(key, p, klen);
    key[klen] = '\0';
    trim(key);

    /* value: strip [units] and trailing comments */
    const char *vstart = eq + 1;
    char buf[MAX_VAL];
    strncpy(buf, vstart, MAX_VAL - 1);
    buf[MAX_VAL - 1] = '\0';

    /* remove trailing comment (# or //) */
    char *cmt = strstr(buf, "//");
    if (cmt) *cmt = '\0';
    cmt = strchr(buf, '#');
    if (cmt) *cmt = '\0';

    /* remove [units] */
    char *bracket = strchr(buf, '[');
    if (bracket) *bracket = '\0';

    strncpy(val, buf, MAX_VAL - 1);
    val[MAX_VAL - 1] = '\0';
    trim(val);
    return 1;
}

static int parse_double_array(const char *s, double *arr, int max_n)
{
    int n = 0;
    char buf[MAX_VAL];
    strncpy(buf, s, MAX_VAL - 1);
    buf[MAX_VAL - 1] = '\0';
    char *tok = strtok(buf, ", \t");
    while (tok && n < max_n) {
        arr[n++] = atof(tok);
        tok = strtok(NULL, ", \t");
    }
    return n;
}

/* ---- parsed v2 config ---- */

typedef struct {
    char forcing_file[MAX_VAL];
    int  verbosity;
    int  num_timesteps;
    double soil_depth_m;
    double soil_b;
    double satdk_m_per_s;
    double satpsi_m;
    double slop;
    double smcmax;
    double wltsmc;
    int    has_wltsmc;
    double alpha_fc;
    double soil_storage_theta;
    double max_gw_storage;
    double cgw;
    double gw_expon;
    double gw_storage;
    double k_nash_subsurface;
    double k_lf;
    double nash_storage_subsurface[MAX_NASH];
    int    num_nash_subsurface;
    char   partitioning_scheme[MAX_VAL];
    double a_xinanjiang;
    double b_xinanjiang;
    double x_xinanjiang;
    double urban_fraction;
    double giuh_ordinates[MAX_GIUH];
    int    num_giuh;
    char   surface_runoff_scheme[MAX_VAL];
    int    has_giuh;
    /* track presence */
    int has_forcing_file;
    int has_soil_depth;
    int has_satdk;
    int has_satpsi;
    int has_smcmax;
    int has_soil_storage;
    int has_partitioning;
} V2Config;

static int read_v2_config(const char *path, V2Config *v2)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "ERROR: Cannot open %s\n", path); return -1; }

    memset(v2, 0, sizeof(*v2));
    v2->verbosity = 0;

    char line[MAX_LINE], key[MAX_KEY], val[MAX_VAL];
    int errors = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (!parse_line(line, key, val)) continue;

        /* check for version tag — already v3 */
        if (ci_eq(key, "cfe_config_version")) {
            fprintf(stderr, "INFO: Config already has cfe_config_version=%s — no migration needed.\n", val);
            fclose(fp);
            return 1; /* signal: already v3 */
        }

        /* map v2 keys */
        if (ci_eq(key, "forcing_file")) {
            strncpy(v2->forcing_file, val, MAX_VAL - 1);
            v2->has_forcing_file = 1;
        }
        else if (ci_eq(key, "verbosity"))              v2->verbosity = atoi(val);
        else if (ci_eq(key, "num_timesteps"))           v2->num_timesteps = atoi(val);
        else if (ci_eq(key, "soil_params.depth") || ci_eq(key, "soil_params.D")) {
            v2->soil_depth_m = atof(val); v2->has_soil_depth = 1;
        }
        else if (ci_eq(key, "soil_params.bb") || ci_eq(key, "soil_params.b"))
            v2->soil_b = atof(val);
        else if (ci_eq(key, "soil_params.satdk")) {
            v2->satdk_m_per_s = atof(val); v2->has_satdk = 1;
        }
        else if (ci_eq(key, "soil_params.satpsi")) {
            v2->satpsi_m = atof(val); v2->has_satpsi = 1;
        }
        else if (ci_eq(key, "soil_params.slop"))        v2->slop = atof(val);
        else if (ci_eq(key, "soil_params.smcmax")) {
            v2->smcmax = atof(val); v2->has_smcmax = 1;
        }
        else if (ci_eq(key, "soil_params.wltsmc")) {
            v2->wltsmc = atof(val); v2->has_wltsmc = 1;
        }
        else if (ci_eq(key, "alpha_fc") || ci_eq(key, "soil_params.alpha_fc"))
            v2->alpha_fc = atof(val);
        else if (ci_eq(key, "soil_storage")) {
            v2->soil_storage_theta = atof(val); v2->has_soil_storage = 1;
        }
        else if (ci_eq(key, "max_gw_storage"))          v2->max_gw_storage = atof(val);
        else if (ci_eq(key, "Cgw"))                     v2->cgw = atof(val);
        else if (ci_eq(key, "expon"))                   v2->gw_expon = atof(val);
        else if (ci_eq(key, "gw_storage"))              v2->gw_storage = atof(val);
        else if (ci_eq(key, "K_nash_subsurface") || ci_eq(key, "K_nash"))
            v2->k_nash_subsurface = atof(val);
        else if (ci_eq(key, "K_lf"))                    v2->k_lf = atof(val);
        else if (ci_eq(key, "nash_storage_subsurface") || ci_eq(key, "nash_storage"))
            v2->num_nash_subsurface = parse_double_array(val, v2->nash_storage_subsurface, MAX_NASH);
        else if (ci_eq(key, "surface_water_partitioning_scheme")) {
            strncpy(v2->partitioning_scheme, val, MAX_VAL - 1);
            v2->has_partitioning = 1;
        }
        else if (ci_eq(key, "a_Xinanjiang_inflection_point_parameter"))
            v2->a_xinanjiang = atof(val);
        else if (ci_eq(key, "b_Xinanjiang_shape_parameter"))
            v2->b_xinanjiang = atof(val);
        else if (ci_eq(key, "x_Xinanjiang_shape_parameter"))
            v2->x_xinanjiang = atof(val);
        else if (ci_eq(key, "urban_decimal_fraction"))
            v2->urban_fraction = atof(val);
        else if (ci_eq(key, "giuh_ordinates")) {
            v2->num_giuh = parse_double_array(val, v2->giuh_ordinates, MAX_GIUH);
            v2->has_giuh = 1;
        }
        else if (ci_eq(key, "surface_runoff_scheme"))
            strncpy(v2->surface_runoff_scheme, val, MAX_VAL - 1);
        /* deprecated keys — silently skip */
        else if (ci_eq(key, "soil_params.expon") ||
                 ci_eq(key, "soil_params.expon_secondary") ||
                 ci_eq(key, "refkdt") ||
                 ci_eq(key, "debug") ||
                 ci_eq(key, "nsubsteps_nash_surface") ||
                 ci_eq(key, "N_nash_surface") ||
                 ci_eq(key, "K_nash_surface") ||
                 ci_eq(key, "nash_storage_surface") ||
                 ci_eq(key, "Kinf_nash_surface") ||
                 ci_eq(key, "retention_depth_nash_surface")) {
            fprintf(stderr, "NOTE: Skipping deprecated v2 key: %s\n", key);
        }
        else {
            fprintf(stderr, "WARNING: Unknown v2 key: %s\n", key);
            errors++;
        }
    }

    fclose(fp);
    if (errors > 0) {
        fprintf(stderr, "ERROR: %d unknown key(s) in v2 config.\n", errors);
        return -1;
    }
    return 0;
}

static int validate_and_convert(V2Config *v2)
{
    /* Check for incompatible surface routing */
    if (strlen(v2->surface_runoff_scheme) > 0 &&
        !ci_eq(v2->surface_runoff_scheme, "GIUH")) {
        if (!v2->has_giuh) {
            fprintf(stderr,
                "ERROR: This config uses surface_runoff_scheme=%s without GIUH ordinates.\n"
                "       Nash Cascade surface routing was removed in CFE v3.\n"
                "       Add giuh_ordinates to the v2 config before migrating, or obtain\n"
                "       GIUH ordinates for this catchment from hydrofabric data.\n",
                v2->surface_runoff_scheme);
            return -1;
        }
        fprintf(stderr,
            "WARNING: surface_runoff_scheme=%s will be replaced with GIUH in v3.\n"
            "         Using provided giuh_ordinates for surface routing.\n",
            v2->surface_runoff_scheme);
    }

    /* Must have GIUH ordinates for v3 */
    if (!v2->has_giuh) {
        fprintf(stderr,
            "ERROR: No giuh_ordinates found. CFE v3 requires GIUH ordinates.\n"
            "       Add giuh_ordinates to the v2 config before migrating.\n");
        return -1;
    }

    /* Required fields */
    if (!v2->has_soil_depth) {
        fprintf(stderr, "ERROR: Missing soil_params.depth\n"); return -1;
    }
    if (!v2->has_satdk) {
        fprintf(stderr, "ERROR: Missing soil_params.satdk\n"); return -1;
    }
    if (!v2->has_satpsi) {
        fprintf(stderr, "ERROR: Missing soil_params.satpsi\n"); return -1;
    }
    if (!v2->has_smcmax) {
        fprintf(stderr, "ERROR: Missing soil_params.smcmax\n"); return -1;
    }

    /* Unit conversions */
    /* satdk: m/s -> cm/h  (× 100 × 3600) */
    v2->satdk_m_per_s = v2->satdk_m_per_s * 100.0 * 3600.0;  /* now cm/h, reuse field */

    /* satpsi: m -> cm  (× 100) */
    v2->satpsi_m = v2->satpsi_m * 100.0;  /* now cm, reuse field */

    /* soil_storage: theta * depth * porosity -> meters */
    if (v2->has_soil_storage) {
        v2->soil_storage_theta = v2->soil_storage_theta * v2->soil_depth_m * v2->smcmax;
        /* now absolute storage in meters, reuse field */
    }

    /* Auto-calculate wilting point if not provided */
    if (!v2->has_wltsmc) {
        v2->wltsmc = wilting_point_from_clapp_hornberger(v2->smcmax, v2->satpsi_m, v2->soil_b);
        if (v2->wltsmc > 0.0) {
            fprintf(stderr, "NOTE: Auto-calculated wilting point = %.6f from Clapp-Hornberger at 15 atm.\n",
                    v2->wltsmc);
        } else {
            v2->wltsmc = 0.066;
            fprintf(stderr, "WARNING: Could not calculate wilting point — using default 0.066.\n");
        }
    }

    /* Default subsurface Nash storage if not provided */
    if (v2->num_nash_subsurface < 1) {
        v2->nash_storage_subsurface[0] = 0.0;
        v2->nash_storage_subsurface[1] = 0.0;
        v2->num_nash_subsurface = 2;
    }

    /* Default partitioning scheme */
    if (!v2->has_partitioning) {
        strncpy(v2->partitioning_scheme, "Schaake", MAX_VAL - 1);
    }

    return 0;
}

static int write_v3_config(const char *path, const V2Config *v2)
{
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "ERROR: Cannot open output %s\n", path); return -1; }

    fprintf(fp, "# CFE v3 Configuration File\n");
    fprintf(fp, "# Migrated from v2 legacy format by cfe_migrate_config\n\n");

    fprintf(fp, "cfe_config_version=3.0[]\n\n");

    fprintf(fp, "# Model Controls\n");
    fprintf(fp, "control_model_timestep_h=1.0[h]\n");
    fprintf(fp, "control_input_forcing_filename=%s\n", v2->forcing_file);
    fprintf(fp, "control_total_num_simulation_timesteps=%d[]\n", v2->num_timesteps);
    fprintf(fp, "control_verbosity=%d[]\n", v2->verbosity);
    fprintf(fp, "control_soil_simulate_freeze_thaw_true_false=FALSE\n");
    fprintf(fp, "control_soil_simulate_discrete_soil_moisture_true_false=FALSE\n");
    fprintf(fp, "control_ET_deepest_root_zone_discretization=4\n");
    fprintf(fp, "control_soil_use_lookup_table_num_points=0\n\n");

    fprintf(fp, "# Catchment Characteristics\n");
    fprintf(fp, "catchment_impervious_fraction_0-1=" FMT "[]\n\n", v2->urban_fraction);

    fprintf(fp, "# Soil Parameters\n");
    fprintf(fp, "soil_depth_m=" FMT "[m]\n", v2->soil_depth_m);
    fprintf(fp, "soil_Clapp_Hornberger_exponent_b=" FMT "[]\n", v2->soil_b);
    fprintf(fp, "soil_sat_hydraulic_conductivity_cm_per_h=" FMT "[cm h-1]\n", v2->satdk_m_per_s);
    fprintf(fp, "soil_sat_capillary_head_cm=" FMT "[cm]\n", v2->satpsi_m);
    fprintf(fp, "soil_effective_porosity=" FMT "[]\n", v2->smcmax);
    // soil_wilting_point_moisture_content removed from v3 — auto-calculated internally
    fprintf(fp, "soil_field_capacity_Pcap_over_Patm_0_1=" FMT "[]\n", v2->alpha_fc);
    fprintf(fp, "soil_to_gw_percolation_rate_limiter_0_to_1=" FMT "[]\n", v2->slop);
    fprintf(fp, "soil_reservoir_rate_const_to_subsurface_lateral_flow=" FMT "[h-1]\n", v2->k_lf);
    fprintf(fp, "soil_ice_content_impervious_threshold=0.0[]\n");

    if (v2->has_soil_storage) {
        fprintf(fp, "state_soil_reservoir_init_storage_m=" FMT "[m]\n\n", v2->soil_storage_theta);
    } else {
        fprintf(fp, "# WARNING: No soil_storage in v2 config — set manually\n");
        fprintf(fp, "state_soil_reservoir_init_storage_m=0.0[m]\n\n");
    }

    fprintf(fp, "# Groundwater Parameters\n");
    fprintf(fp, "gw_reservoir_max_storage_m=" FMT "[m]\n", v2->max_gw_storage);
    fprintf(fp, "gw_discharge_coeff_m_per_timestep=" FMT "[m h-1]\n", v2->cgw);
    fprintf(fp, "gw_discharge_exponent=" FMT "[]\n", v2->gw_expon);
    fprintf(fp, "state_gw_reservoir_init_storage_m=" FMT "[m]\n\n", v2->gw_storage);

    fprintf(fp, "# Partitioning Scheme\n");
    fprintf(fp, "partitioning_scheme_name=%s\n", v2->partitioning_scheme);
    fprintf(fp, "partitioning_Xinanjiang_tension_water_inflection_point=" FMT "[]\n", v2->a_xinanjiang);
    fprintf(fp, "partitioning_Xinanjiang_tension_water_soil_moist_distrib_exponent=" FMT "[]\n", v2->b_xinanjiang);
    fprintf(fp, "partitioning_Xinanjiang_free_water_soil_moist_distrib_exponent=" FMT "[]\n\n", v2->x_xinanjiang);

    fprintf(fp, "# Surface Routing (GIUH only in CFE v3)\n");
    fprintf(fp, "surface_routing_num_giuh_ordinates=%d\n", v2->num_giuh);
    fprintf(fp, "surface_routing_giuh_ordinates=");
    for (int i = 0; i < v2->num_giuh; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, FMT, v2->giuh_ordinates[i]);
    }
    fprintf(fp, "[]\n");

    /* GIUH convolution queue: zeros matching ordinate count */
    fprintf(fp, "state_surface_routing_init_giuh_convolution_queue_m=");
    for (int i = 0; i < v2->num_giuh; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, "0.0");
    }
    fprintf(fp, "[m]\n\n");

    fprintf(fp, "# Subsurface Lateral Flow (Nash Cascade)\n");
    fprintf(fp, "subsurface_routing_nash_reservoir_time_constant_k=" FMT "[h-1]\n", v2->k_nash_subsurface);
    fprintf(fp, "state_subsurface_routing_init_nash_cascade_storage_m=");
    for (int i = 0; i < v2->num_nash_subsurface; i++) {
        if (i > 0) fprintf(fp, ", ");
        fprintf(fp, FMT, v2->nash_storage_subsurface[i]);
    }
    fprintf(fp, "[m]\n");

    fclose(fp);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Usage: cfe_migrate_config <input_v2_config> <output_v3_config>\n");
        return 1;
    }

    V2Config v2;
    int rc = read_v2_config(argv[1], &v2);
    if (rc == 1) return 0;  /* already v3 */
    if (rc != 0) return 1;

    if (validate_and_convert(&v2) != 0) return 1;
    if (write_v3_config(argv[2], &v2) != 0) return 1;

    fprintf(stderr, "Migrated: %s -> %s\n", argv[1], argv[2]);
    return 0;
}
