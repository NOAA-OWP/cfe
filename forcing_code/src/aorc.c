#ifndef AORC_C
#define AORC_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//local includes
#include "../include/aorc.h"
#include "../include/aorc_tools.h"

extern void alloc_aorc_model(aorc_model *model) {
    // TODO: *******************
}

extern void free_aorc_model(aorc_model *model) {
    // TODO: ******************
}

// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
extern int run_aorc(aorc_model* model)
{

  // populate the evapotranspiration forcing data structure:
  //---------------------------------------------------------------------------------------------------------------
  model->other_forcing.air_temperature_C = model->forcing_data_air_temperature_2m_K[model->bmi.current_step] - TK;//convert to C
  model->other_forcing.relative_humidity_percent     = (double)-99.9; // this negative number means use specific humidity
  model->other_forcing.specific_humidity_2m_kg_per_kg = model->forcing_data_precip_kg_per_m2[model->bmi.current_step];
  model->other_forcing.air_pressure_Pa    = model->forcing_data_surface_pressure_Pa[model->bmi.current_step];
  model->other_forcing.wind_speed_m_per_s = hypot(model->forcing_data_u_wind_speed_10m_m_per_s[model->bmi.current_step],
                                         model->forcing_data_v_wind_speed_10m_m_per_s[model->bmi.current_step]);                 

  model->aorc.precip_kg_per_m2               =  model->forcing_data_precip_kg_per_m2[model->bmi.current_step];
  model->aorc.incoming_longwave_W_per_m2     =  model->forcing_data_incoming_longwave_W_per_m2[model->bmi.current_step];
  model->aorc.incoming_shortwave_W_per_m2    =  model->forcing_data_incoming_shortwave_W_per_m2[model->bmi.current_step];
  model->aorc.surface_pressure_Pa            =  model->forcing_data_surface_pressure_Pa[model->bmi.current_step];
  model->aorc.specific_humidity_2m_kg_per_kg =  model->forcing_data_specific_humidity_2m_kg_per_kg[model->bmi.current_step];
  // jframe: not sure if this belongs here or not, but it needs to happen somewhere.
  model->other_forcing.specific_humidity_2m_kg_per_kg =  model->aorc.specific_humidity_2m_kg_per_kg;

  model->aorc.air_temperature_2m_K           =  model->forcing_data_air_temperature_2m_K[model->bmi.current_step];
  model->aorc.u_wind_speed_10m_m_per_s       =  model->forcing_data_u_wind_speed_10m_m_per_s[model->bmi.current_step];
  model->aorc.v_wind_speed_10m_m_per_s       =  model->forcing_data_v_wind_speed_10m_m_per_s[model->bmi.current_step];
  model->aorc.latitude                       =  model->solar_params.latitude_degrees;
  model->aorc.longitude                      =  model->solar_params.longitude_degrees;

  // wind speed was measured at 10.0 m height, so we need to calculate the wind speed at 2.0m
  double numerator=log(2.0/model->aorc_params.zero_plane_displacement_height_m);
  double denominator=log(model->aorc_params.wind_speed_measurement_height_m/model->aorc_params.zero_plane_displacement_height_m);
  model->other_forcing.wind_speed_m_per_s = model->other_forcing.wind_speed_m_per_s*numerator/denominator;  // this is the 2 m value
  model->aorc_params.wind_speed_measurement_height_m=2.0;  // change because we converted from 10m to 2m height.
  // transfer aorc forcing data into our data structure for surface radiation calculations
  model->surf_rad_forcing.incoming_shortwave_radiation_W_per_sq_m = (double)model->aorc.incoming_shortwave_W_per_m2;
  model->surf_rad_forcing.incoming_longwave_radiation_W_per_sq_m  = (double)model->aorc.incoming_longwave_W_per_m2; 
  model->surf_rad_forcing.air_temperature_C                       = (double)model->aorc.air_temperature_2m_K-TK;

  // compute relative humidity from specific humidity..
  double saturation_vapor_pressure_Pa = aorc_calc_air_saturation_vapor_pressure_Pa(model->surf_rad_forcing.air_temperature_C);
  double actual_vapor_pressure_Pa = (double)model->aorc.specific_humidity_2m_kg_per_kg*(double)model->aorc.surface_pressure_Pa/0.622;

  model->surf_rad_forcing.relative_humidity_percent = 100.0*actual_vapor_pressure_Pa/saturation_vapor_pressure_Pa;
  // sanity check the resulting value.  Should be less than 100%.  Sometimes air can be supersaturated.
  if(100.0< model->surf_rad_forcing.relative_humidity_percent) model->surf_rad_forcing.relative_humidity_percent = 99.0;


  aorc_calculate_solar_radiation(model);
  
  if (model->bmi.verbose > 1)
    printf("calculate the net radiation\n");
    // NOTE don't call this function use_aerodynamic_method option is TRUE
  model->other_forcing.net_radiation_W_per_sq_m=aorc_calculate_net_radiation_W_per_sq_m(model);

  return 0;
}

//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
void aorc_setup(aorc_model* model)
{
  model->other_forcing.canopy_resistance_sec_per_m   = 50.0; // TODO: from plant growth model
  model->other_forcing.water_temperature_C           = 15.5; // TODO: from soil or lake thermal model
  model->other_forcing.ground_heat_flux_W_per_sq_m   = -10.0;    // TODO from soil thermal model.  Negative denotes downward.

  // these values are needed if we don't have incoming longwave radiation measurements.
  model->surf_rad_forcing.incoming_shortwave_radiation_W_per_sq_m     = 440.1;     // must come from somewhere
  model->surf_rad_forcing.incoming_longwave_radiation_W_per_sq_m      = -1.0e+05;  // this huge negative value tells to calc.
  model->surf_rad_forcing.air_temperature_C                           = 15.0;      // from some forcing data file
  model->surf_rad_forcing.relative_humidity_percent                   = 63.0;      // from some forcing data file
  model->surf_rad_forcing.ambient_temperature_lapse_rate_deg_C_per_km = 6.49;      // ICAO standard atmosphere lapse rate
  model->surf_rad_forcing.cloud_cover_fraction                        = 0.6;       // from some forcing data file
  model->surf_rad_forcing.cloud_base_height_m                         = 2500.0/3.281; // assumed 2500 ft.

  // these values are needed if we don't have incoming longwave radiation measurements.
  model->surf_rad_forcing.incoming_shortwave_radiation_W_per_sq_m     = 440.1;     // must come from somewhere
  model->surf_rad_forcing.incoming_longwave_radiation_W_per_sq_m      = -1.0e+05;  // this huge negative value tells to calc.
  model->surf_rad_forcing.air_temperature_C                           = 15.0;      // from some forcing data file
  model->surf_rad_forcing.relative_humidity_percent                   = 63.0;      // from some forcing data file
  model->surf_rad_forcing.ambient_temperature_lapse_rate_deg_C_per_km = 6.49;      // ICAO standard atmosphere lapse rate
  model->surf_rad_forcing.cloud_cover_fraction                        = 0.6;       // from some forcing data file
  model->surf_rad_forcing.cloud_base_height_m                         = 2500.0/3.281; // assumed 2500 ft.
    
  // Surface radiation forcing parameter values that must come from other models
  //---------------------------------------------------------------------------------------------------------------
  model->surf_rad_forcing.surface_skin_temperature_C = 12.0;  // TODO from soil thermal model or vegetation model.

  model->solar_options.cloud_base_height_known=0;  // set to TRUE if the solar_forcing.cloud_base_height_m is known.

  // ### FORCING ###
  model->surf_rad_forcing.cloud_cover_fraction         =   0.5;   // THESE VALUES ARE FOR THE UNIT TEST 
  model->surf_rad_forcing.atmospheric_turbidity_factor =   2.0;   // 2.0 = clear mountain air, 5.0= smoggy air
  model->surf_rad_forcing.day_of_year                  =  208;    // THESE VALUES ARE FOR THE UNIT TEST
  model->surf_rad_forcing.zulu_time                  =  20.567; // THESE VALUES ARE FOR THE UNIT TEST

  return;
}


void aorc_unit_tests(aorc_model* model)
{
  printf("\n #----------- BEGIN UNIT TESTS   ---------------# \n");
  
  printf("\n #-----------       UNIT TEST    ---------------# \n");
  printf("shortwave radiation is %lf W s-1,\n and should be: 117.1 W -1 \n",
      model->aorc.incoming_shortwave_W_per_m2);

  printf("\n #------------  END UNIT TESTS   ---------------# \n");
  printf("\n");
  printf("\n");
  return;
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* ALL THE STUFF BELOW HERE IS JUST UTILITY MEMORY AND TIME FUNCTION CODE */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/


/*####################################################################*/
/*########################### PARSE LINE #############################*/
/*####################################################################*/
void parse_aorc_line_aorc(char *theString,long *year,long *month, long *day,long *hour,long *minute, double *second,
                struct aorc_forcing_data *aorc)
{
    int start,end;
    int wordlen;
    char theWord[150];
    
    //len=strlen(theString);
    
    start=0; /* begin at the beginning of theString */
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    *year=atol(theWord);
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    *month=atol(theWord);
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    *day=atol(theWord);
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    *hour=atol(theWord);
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    *minute=atol(theWord);
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    *second=(double)atof(theWord);
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->precip_kg_per_m2=atof(theWord);
    //printf("%s, %s, %lf, %lf\n", theString, theWord, (double)atof(theWord), aorc->precip_kg_per_m2);
                  
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->incoming_longwave_W_per_m2=(double)atof(theWord);   
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->incoming_shortwave_W_per_m2=(double)atof(theWord);   
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->surface_pressure_Pa=(double)atof(theWord);           
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->specific_humidity_2m_kg_per_kg=(double)atof(theWord);
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->air_temperature_2m_K=(double)atof(theWord);          
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->u_wind_speed_10m_m_per_s=(double)atof(theWord);      
    
    get_word_aorc(theString,&start,&end,theWord,&wordlen);
    aorc->v_wind_speed_10m_m_per_s=(double)atof(theWord);      

    //----------------------------------------------    
    // Is this needed?  It has not been set. (SDP)
    //----------------------------------------------
    aorc->time = -9999.0;
    //######################################################
     
    return;
    }

/*####################################################################*/
/*############################## GET WORD ############################*/
/*####################################################################*/
void get_word_aorc(char *theString, int *start, int *end, char *theWord, int *wordlen) {
    int i, lenny, j;
    lenny = strlen(theString);

    while (theString[*start] == ' ' || theString[*start] == '\t') {
        (*start)++;
    };

    j = 0;
    for (i = (*start); i < lenny; i++) {
        if (theString[i] != ' ' && theString[i] != '\t' && theString[i] != ',' && theString[i] != ':' &&
            theString[i] != '/') {
            theWord[j] = theString[i];
            j++;
        }
        else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
             && (theString[i]=='-' && theString[i-1]==',') )
        {
            theWord[j]=theString[i];
            j++;
        }
        else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
             && (theString[i]=='-' && theString[i-1]=='e') )
        {
            theWord[j]=theString[i];
            j++;
        }
        else {
            break;
        }
    }
    theWord[j] = '\0';
    *start = i + 1;
    *wordlen = strlen(theWord);
    return;
}

/****************************************/
void itwo_alloc_aorc(int ***array, int rows, int cols) {
    int i, frows, fcols, numgood = 0;
    int error = 0;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
        exit(-9);
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (int **) malloc(frows * sizeof(int *));
    if (*array) {
        memset((*array), 0, frows * sizeof(int *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (int *) malloc(fcols * sizeof(int));
            if ((*array)[i] == NULL) {
                error = 1;
                numgood = i;
                i = frows;
            } else
                memset((*array)[i], 0, fcols * sizeof(int));
        }
    }
    return;
}


void dtwo_alloc_aorc(double ***array, int rows, int cols) {
    int i, frows, fcols, numgood = 0;
    int error = 0;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
        exit(-9);
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (double **) malloc(frows * sizeof(double *));
    if (*array) {
        memset((*array), 0, frows * sizeof(double *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (double *) malloc(fcols * sizeof(double));
            if ((*array)[i] == NULL) {
                error = 1;
                numgood = i;
                i = frows;
            } else
                memset((*array)[i], 0, fcols * sizeof(double));
        }
    }
    return;
}


void d_alloc_aorc(double **var, int size) {
    size++;  /* just for safety */

    *var = (double *) malloc(size * sizeof(double));
    if (*var == NULL) {
        printf("Problem allocating memory for array in d_alloc.\n");
        return;
    } else
        memset(*var, 0, size * sizeof(double));
    return;
}

void i_alloc_aorc(int **var, int size) {
    size++;  /* just for safety */

    *var = (int *) malloc(size * sizeof(int));
    if (*var == NULL) {
        printf("Problem allocating memory in i_alloc\n");
        return;
    } else
        memset(*var, 0, size * sizeof(int));
    return;
}

/*
 * convert Gregorian days to Julian date
 *
 * Modify as needed for your application.
 *
 * The Julian day starts at noon of the Gregorian day and extends
 * to noon the next Gregorian day.
 *
 */
/*
** Takes a date, and returns a Julian day. A Julian day is the number of
** days since some base date  (in the very distant past).
** Handy for getting date of x number of days after a given Julian date
** (use jdate to get that from the Gregorian date).
** Author: Robert G. Tantzen, translator: Nat Howard
** Translated from the algol original in Collected Algorithms of CACM
** (This and jdate are algorithm 199).
*/


double greg_2_jul_aorc(
        long year,
        long mon,
        long day,
        long h,
        long mi,
        double se) {
    long m = mon, d = day, y = year;
    long c, ya, j;
    double seconds = h * 3600.0 + mi * 60 + se;

    if (m > 2)
        m -= 3;
    else {
        m += 9;
        --y;
    }
    c = y / 100L;
    ya = y - (100L * c);
    j = (146097L * c) / 4L + (1461L * ya) / 4L + (153L * m + 2L) / 5L + d + 1721119L;
    if (seconds < 12 * 3600.0) {
        j--;
        seconds += 12.0 * 3600.0;
    } else {
        seconds = seconds - 12.0 * 3600.0;
    }
    return (j + (seconds / 3600.0) / 24.0);
}

/* Julian date converter. Takes a julian date (the number of days since
** some distant epoch or other), and returns an int pointer to static space.
** ip[0] = month;
** ip[1] = day of month;
** ip[2] = year (actual year, like 1977, not 77 unless it was  77 a.d.);
** ip[3] = day of week (0->Sunday to 6->Saturday)
** These are Gregorian.
** Copied from Algorithm 199 in Collected algorithms of the CACM
** Author: Robert G. Tantzen, Translator: Nat Howard
**
** Modified by FLO 4/99 to account for nagging round off error 
**
*/
void calc_date_aorc(double jd, long *y, long *m, long *d, long *h, long *mi,
               double *sec) {
    static int ret[4];

    long j;
    double tmp;
    double frac;

    j = (long) jd;
    frac = jd - j;

    if (frac >= 0.5) {
        frac = frac - 0.5;
        j++;
    } else {
        frac = frac + 0.5;
    }

    ret[3] = (j + 1L) % 7L;
    j -= 1721119L;
    *y = (4L * j - 1L) / 146097L;
    j = 4L * j - 1L - 146097L * *y;
    *d = j / 4L;
    j = (4L * *d + 3L) / 1461L;
    *d = 4L * *d + 3L - 1461L * j;
    *d = (*d + 4L) / 4L;
    *m = (5L * *d - 3L) / 153L;
    *d = 5L * *d - 3 - 153L * *m;
    *d = (*d + 5L) / 5L;
    *y = 100L * *y + j;
    if (*m < 10)
        *m += 3;
    else {
        *m -= 9;
        *y = *y + 1; /* Invalid use: *y++. Modified by Tony */
    }

    /* if (*m < 3) *y++; */
    /* incorrectly repeated the above if-else statement. Deleted by Tony.*/

    tmp = 3600.0 * (frac * 24.0);
    *h = (long) (tmp / 3600.0);
    tmp = tmp - *h * 3600.0;
    *mi = (long) (tmp / 60.0);
    *sec = tmp - *mi * 60.0;
}





#endif  // AORC_C
