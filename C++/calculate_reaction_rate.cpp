#ifndef CALC_RATE_H
#define CALC_RATE_H

#include "math.h"
#include "types.hpp"
#include "constants.hpp"
#include "cvodes/cvodes.h"

namespace CALC_RATE {


inline TYPES::DTP_FLOAT thermal_velocity_CGS(const TYPES::DTP_FLOAT T_CGS, const TYPES::DTP_FLOAT massnum) {
  return sqrt((8.0/CONST::PI)
              * CONST::phy_kBoltzmann_CGS * T_CGS
              / (massnum * CONST::phy_mProton_CGS));
}


TYPES::DTP_FLOAT rate_adsorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * r.abc[0] // Act as sticking efficiency
       * p.get_n_gas()
       * p.get_dust2gas_num()
       * p.get_dust_crosssec()
       * thermal_velocity_CGS(p.get_T_gas(),
                              s.massSpecies.at(r.idxReactants[0]));
}


std::vector<TYPES::DTP_FLOAT> drdy_adsorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(r.abc[0]
            * p.get_n_gas()
            * p.get_dust2gas_num()
            * p.get_dust_crosssec()
            * thermal_velocity_CGS(p.get_T_gas(),
                                   s.massSpecies.at(r.idxReactants[0])));
  return drdy;
}


TYPES::DTP_FLOAT rate_desorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * s.vibFreqs.at(r.idxReactants[0])
       * r.abc[0]
       * (exp(-r.abc[2]/p.get_T_dust())
        + p.get_chi_cosmicray()
        * CONST::phy_cosmicray_desorption_factor
        * exp(-r.abc[2] / CONST::phy_cosmicray_desorption_T));
}


std::vector<TYPES::DTP_FLOAT> drdy_desorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(
      s.vibFreqs.at(r.idxReactants[0])
    * r.abc[0]
    * (exp(-r.abc[2]/p.get_T_dust())
     + p.get_chi_cosmicray()
     * CONST::phy_cosmicray_desorption_factor
     * exp(-r.abc[2] / CONST::phy_cosmicray_desorption_T)));
  return drdy;
}


TYPES::DTP_FLOAT rate_ion_neutral(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * NV_Ith_S(y, r.idxReactants[1])
       * r.abc[0]
       * pow(p.get_T_gas()/3e2, r.abc[1])
       * exp(-r.abc[2] / p.get_T_gas());
}


std::vector<TYPES::DTP_FLOAT> drdy_ion_neutral(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  TYPES::DTP_FLOAT tmp = r.abc[0] * pow(p.get_T_gas()/3e2, r.abc[1])
                       * exp(-r.abc[2] / p.get_T_gas());
  drdy.push_back(tmp * NV_Ith_S(y, r.idxReactants[1]));
  drdy.push_back(tmp * NV_Ith_S(y, r.idxReactants[0]));
  return drdy;
}


TYPES::DTP_FLOAT rate_cosmicray_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * r.abc[0] * p.get_chi_cosmicray();
}


std::vector<TYPES::DTP_FLOAT> drdy_cosmicray_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(r.abc[0] * p.get_chi_cosmicray());
  return drdy;
}


TYPES::DTP_FLOAT rate_cosmicray_induced_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * r.abc[0] * p.get_chi_cosmicray()
       * pow(p.get_T_gas()/3e2, r.abc[1])
       * r.abc[2] / (1.0 - p.get_dust_albedo());
}


std::vector<TYPES::DTP_FLOAT> drdy_cosmicray_induced_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(r.abc[0] * p.get_chi_cosmicray()
               * pow(p.get_T_gas()/3e2, r.abc[1])
               * r.abc[2] / (1.0 - p.get_dust_albedo()));
  return drdy;
}


TYPES::DTP_FLOAT rate_photoionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * r.abc[0] * p.get_G0_UV() * exp(-r.abc[2] * p.get_Av());
}


std::vector<TYPES::DTP_FLOAT> drdy_photoionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(r.abc[0] * p.get_G0_UV()
               * exp(-r.abc[2] * p.get_Av()));
  return drdy;
}


int assignReactionHandlers(TYPES::User_data& user_data) {
  (*user_data.rate_calculators)[1]  = rate_cosmicray_ionization;
  (*user_data.drdy_calculators)[1]  = drdy_cosmicray_ionization;
  (*user_data.rate_calculators)[71] = rate_cosmicray_ionization;
  (*user_data.drdy_calculators)[71] = drdy_cosmicray_ionization;
  (*user_data.rate_calculators)[2]  = rate_cosmicray_induced_ionization;
  (*user_data.drdy_calculators)[2]  = drdy_cosmicray_induced_ionization;
  (*user_data.rate_calculators)[72] = rate_cosmicray_induced_ionization;
  (*user_data.drdy_calculators)[72] = drdy_cosmicray_induced_ionization;
  (*user_data.rate_calculators)[3]  = rate_photoionization;
  (*user_data.drdy_calculators)[3]  = drdy_photoionization;
  (*user_data.rate_calculators)[73] = rate_photoionization;
  (*user_data.drdy_calculators)[73] = drdy_photoionization;
  (*user_data.rate_calculators)[5]  = rate_ion_neutral;
  (*user_data.drdy_calculators)[5]  = drdy_ion_neutral;
  (*user_data.rate_calculators)[53] = rate_ion_neutral;
  (*user_data.drdy_calculators)[53] = drdy_ion_neutral;
  (*user_data.rate_calculators)[61] = rate_adsorption;
  (*user_data.drdy_calculators)[61] = drdy_adsorption;
  (*user_data.rate_calculators)[62] = rate_desorption;
  (*user_data.drdy_calculators)[62] = drdy_desorption;
  return 0;
}

}

#endif //CALC_RATE_H
