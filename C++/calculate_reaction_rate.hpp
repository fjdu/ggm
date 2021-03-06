#ifndef CALC_RATE_H
#define CALC_RATE_H

#include "math.h"
#include <algorithm>
#include "types.hpp"
#include "constants.hpp"


namespace CALC_RATE {


extern inline TYPES::DTP_FLOAT thermal_velocity_CGS(
    const TYPES::DTP_FLOAT T_CGS,
    const TYPES::DTP_FLOAT massnum);


extern TYPES::DTP_FLOAT rate_adsorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_desorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_ion_neutral(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_cosmicray_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_cosmicray_induced_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_photoionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern inline TYPES::DTP_FLOAT calc_cross_surf_barrier_prob(
    const TYPES::PhyParams& p, TYPES::Reaction& r);


extern TYPES::DTP_FLOAT rate_surface_AA(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_surface_AB(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_surf2mant(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern void update_surfmant(
    const TYPES::DTP_FLOAT& t,
    double *y,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_mant2surf(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_iongrain(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_photodesorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern TYPES::DTP_FLOAT rate_dummy(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m);


extern int assignReactionHandlers(TYPES::User_data& user_data);

}

#endif //CALC_RATE_H
