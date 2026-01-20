extern "C" {
#include "bmi_serialization.h"
#include "bmi_cfe.h"
#include "logger.h"
}
#include "vecbuf.hpp"

#include <sstream>

#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


class CfeSerializer {
    public:
        CfeSerializer(Bmi* bmi);
        cfe_state_struct* data;
    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int version);
};


CfeSerializer::CfeSerializer(Bmi* bmi) {
    void* data = bmi->data;
    this->data = (cfe_state_struct*)data;
}


template<class Archive>
void CfeSerializer::serialize(Archive& ar, const unsigned int version) {
    using boost::serialization::make_array;
    cfe_state_struct* state = this->data;
    
    // base
    ar & state->current_time_step; // update
    ar & state->timestep_rainfall_input_m; // update
    ar & state->nwm_ponded_depth_m; // update
    ar & state->soil_reservoir_storage_deficit_m; // cfe
    ar & *state->infiltration_excess_m; // cfe; stored as pointer on struct
    ar & state->infiltration_depth_m; // cfe
    ar & *state->flux_perc_m; // cfe; stored as pointer on struct
    ar & *state->flux_lat_m; // cfe; stored as pointer on struct
    ar & state->gw_reservoir_storage_deficit_m; // cfe
    ar & *state->flux_from_deep_gw_to_chan_m; // cfe; stored as pointer on struct
    ar & *state->flux_direct_runoff_m; // cfe; stored as pointer on struct
    ar & *state->flux_nash_lateral_runoff_m; // cfe, nash_cascade_subsurface; stored as pointer on struct
    ar & *state->flux_Qout_m; // cfe; stored as pointer on struct
    ar & make_array(state->nash_storage_subsurface, // nash_cascade_subsurface
                    state->N_nash_subsurface); // config

    // conceptual_reservoir
    ar & state->soil_reservoir.storage_m; // cfe, et_from_soil
    ar & state->soil_reservoir.storage_change_m; // cfe
    ar & make_array(state->soil_reservoir.smc_profile, // et_from_soil
                    state->soil_reservoir.n_soil_layers); // config
    ar & state->soil_reservoir.ice_fraction_schaake; // cfe
    ar & state->soil_reservoir.ice_fraction_xinanjiang; // cfe
    ar & state->gw_reservoir.storage_m; // cfe

    // nash_cascade_parameters
    if (state->nash_surface_params.nash_storage != NULL) { // will be null or malloc'd based on config
        ar & make_array(state->nash_surface_params.nash_storage, // et_from_retention_depth, nash_cascade_surface
                        state->nash_surface_params.N_nash); // config
    }

    // evapotransporation_struct
    ar & state->et_struct.potential_et_m_per_timestep; // cfe
    ar & state->et_struct.reduced_potential_et_m_per_timestep; // cfe, et_from_soil, et_from_rainfall
    ar & state->et_struct.actual_et_from_rain_m_per_timestep; // cfe
    ar & state->et_struct.actual_et_from_soil_m_per_timestep; // cfe, et_from_soil
    ar & state->et_struct.actual_et_m_per_timestep; // cfe
    ar & state->et_struct.potential_et_m_per_s; // cfe

    // massball, all set in cfe
    ar & state->vol_struct.volin; // set in update
    ar & state->vol_struct.vol_et_from_rain;
    ar & state->vol_struct.vol_et_to_atm;
    ar & state->vol_struct.volout;
    ar & state->vol_struct.vol_et_from_retention_depth;
    ar & state->vol_struct.vol_out_surface;
    ar & state->vol_struct.vol_et_from_soil;
    ar & state->vol_struct.vol_runoff;
    ar & state->vol_struct.vol_infilt;
    ar & state->vol_struct.vol_to_soil;
    ar & state->vol_struct.vol_to_gw;
    ar & state->vol_struct.vol_soil_to_gw;
    ar & state->vol_struct.vol_soil_to_lat_flow;
    ar & state->vol_struct.vol_from_gw;
    ar & state->vol_struct.vol_out_surface;
    ar & state->vol_struct.vol_in_nash;
    ar & state->vol_struct.vol_out_nash;

    if (state->surface_runoff_scheme == GIUH) { // config
        ar & make_array(state->runoff_queue_m_per_timestep, // giuh_convolution_integral
                        state->num_giuh_ordinates + 1); // config
    } 
    /// state->nash_surface_params.runon_infiltration not used for input or GetValue
    // else if (state->surface_runoff_scheme == NASH_CASCADE) {
    //     ar & state->nash_surface_params.runon_infiltration;
    // }
}


extern "C" {

int free_serialized_cfe(Bmi* bmi) {
    cfe_state_struct* model = (cfe_state_struct*)bmi->data;
    if (model->serialized != NULL) {
        free(model->serialized);
        model->serialized = NULL;
    }
    model->serialized_length = 0;
    return BMI_SUCCESS;
}

int load_serialized_cfe(Bmi* bmi, const char* data) {
    CfeSerializer serializer(bmi);
    std::istringstream stream(data);
    boost::archive::binary_iarchive archive(stream);
    try {
        archive >> serializer;
        return BMI_SUCCESS;
    } catch (const std::exception &e) {
        Log(LogLevel::SEVERE, "Deserializing CFE encountered an error: %s", e.what());
        return BMI_FAILURE;
    }
}

int new_serialized_cfe(Bmi* bmi) {
    CfeSerializer serializer(bmi);
    vecbuf<char> stream;
    boost::archive::binary_oarchive archive(stream);
    try {
        archive << serializer;
    } catch (const std::exception &e) {
        Log(LogLevel::SEVERE, "Serializing CFE encountered an error: %s", e.what());
        return BMI_FAILURE;
    }
    // copy serialized data into cfe data
    cfe_state_struct* model = (cfe_state_struct*)bmi->data;
    // clear previous data if it exists
    if (model->serialized != NULL) {
        free(model->serialized);
    }
    // set size and allocate memory
    model->serialized_length = stream.size();
    model->serialized = (char*)malloc(model->serialized_length);
    // make sure memory could be allocated
    if (model->serialized == NULL) {
        Log(LogLevel::SEVERE, "CFE failed to allocate memory for serializing.");
        model->serialized_length = 0;
        return BMI_FAILURE;
    }
    // copy stream data to new allocation
    memcpy(model->serialized, stream.data(), model->serialized_length);
    return BMI_SUCCESS;
}

}
