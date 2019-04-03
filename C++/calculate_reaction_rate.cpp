#ifndef CALC_RATE_H
#define CALC_RATE_H

#include "math.h"
#include <algorithm>
#include "types.hpp"
#include "constants.hpp"
#include "cvodes/cvodes.h"

namespace CALC_RATE {


inline TYPES::DTP_FLOAT thermal_velocity_CGS(
    const TYPES::DTP_FLOAT T_CGS,
    const TYPES::DTP_FLOAT massnum) {
  return sqrt((8.0/CONST::PI)
              * CONST::phy_kBoltzmann_CGS * T_CGS
              / (massnum * CONST::phy_mProton_CGS));
}


TYPES::DTP_FLOAT rate_adsorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
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
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(
         r.abc[0]
       * p.get_n_gas()
       * p.get_dust2gas_num()
       * p.get_dust_crosssec()
       * thermal_velocity_CGS(p.get_T_gas(),
                              s.massSpecies.at(r.idxReactants[0])));
  //std::cout << s.idx2name.at(r.idxReactants[0]) << " "
  //          << drdy[0] << std::endl;
  return drdy;
}


TYPES::DTP_FLOAT rate_desorption(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
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
    const TYPES::Species& s,
    TYPES::OtherData& m)
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
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * NV_Ith_S(y, r.idxReactants[1])
       * p.get_n_gas()
       * r.abc[0]
       * pow(p.get_T_gas()/3e2, r.abc[1])
       * exp(-r.abc[2] / p.get_T_gas());
}


std::vector<TYPES::DTP_FLOAT> drdy_ion_neutral(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  TYPES::DTP_FLOAT tmp = p.get_n_gas() * r.abc[0]
                       * pow(p.get_T_gas()/3e2, r.abc[1])
                       * exp(-r.abc[2] / p.get_T_gas());
  drdy.push_back(tmp * NV_Ith_S(y, r.idxReactants[1]));
  drdy.push_back(tmp * NV_Ith_S(y, r.idxReactants[0]));
  //std::cout << s.idx2name.at(r.idxReactants[0]) << " "
  //          << s.idx2name.at(r.idxReactants[1]) << " "
  //          << tmp << std::endl;
  return drdy;
}


TYPES::DTP_FLOAT rate_cosmicray_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * r.abc[0] * p.get_chi_cosmicray();
}


std::vector<TYPES::DTP_FLOAT> drdy_cosmicray_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(r.abc[0] * p.get_chi_cosmicray());
  //std::cout << s.idx2name.at(r.idxReactants[0]) << " "
  //          << drdy[0] << std::endl;
  return drdy;
}


TYPES::DTP_FLOAT rate_cosmicray_induced_ionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
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
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(r.abc[0] * p.get_chi_cosmicray()
               * pow(p.get_T_gas()/3e2, r.abc[1])
               * r.abc[2] / (1.0 - p.get_dust_albedo()));
  //std::cout << s.idx2name.at(r.idxReactants[0]) << " "
  //          << drdy[0] << std::endl;
  return drdy;
}


TYPES::DTP_FLOAT rate_photoionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  return NV_Ith_S(y, r.idxReactants[0])
       * r.abc[0] * p.get_G0_UV() * exp(-r.abc[2] * p.get_Av());
}


std::vector<TYPES::DTP_FLOAT> drdy_photoionization(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(r.abc[0] * p.get_G0_UV()
               * exp(-r.abc[2] * p.get_Av()));
  std::cout << s.idx2name.at(r.idxReactants[0]) << " "
            << drdy[0] << std::endl;
  return drdy;
}


inline TYPES::DTP_FLOAT calc_cross_surf_barrier_prob(
    const TYPES::PhyParams& p, const TYPES::Reaction& r)
{
  TYPES::DTP_FLOAT tmp = 1.0;
  if ((r.abc[2] > 0.0) && (r.Trange[0] > 0.0)) {
    tmp *=
      exp(std::max(-r.abc[2]/p.get_T_dust(),
                   -2.0 * r.abc[1] * CONST::phy_DiffBarrierWidth_CGS
                   / CONST::phy_hbarPlanck_CGS
                   * sqrt(2.0 * r.Trange[0] * CONST::phy_mProton_CGS
                        * CONST::phy_kBoltzmann_CGS * r.abc[2])));
  }
  return tmp;
}


TYPES::DTP_FLOAT rate_surface_AA(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  int i0 = r.idxReactants[0];
  //std::cout << i0 << s.idx2name.at(i0) << std::endl;
  TYPES::DTP_FLOAT tmp =
         std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/p.get_T_dust()),
                  s.quantMobilities.at(i0))
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density())
       / p.get_dust2gas_num() * r.abc[0] * calc_cross_surf_barrier_prob(p, r);
  return NV_Ith_S(y, i0)
       * NV_Ith_S(y, i0)
       * tmp;
}


std::vector<TYPES::DTP_FLOAT> drdy_surface_AA(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  int i0 = r.idxReactants[0];
  TYPES::DTP_FLOAT tmp =
         std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/p.get_T_dust()),
                  s.quantMobilities.at(i0))
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density())
       / p.get_dust2gas_num() * r.abc[0] * calc_cross_surf_barrier_prob(p, r);
  drdy.push_back(tmp * NV_Ith_S(y, r.idxReactants[1]));
  drdy.push_back(tmp * NV_Ith_S(y, r.idxReactants[0]));
  return drdy;
}


TYPES::DTP_FLOAT rate_surface_AB(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  int i0 = r.idxReactants[0], i1 = r.idxReactants[1];
  //std::cout << i0 << s.idx2name.at(i0) << " " << s.idx2name.at(i1) << std::endl;
  TYPES::DTP_FLOAT tmp =
         (std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/p.get_T_dust()),
                   s.quantMobilities.at(i0))
        + std::max(s.vibFreqs.at(i1) * exp(-s.diffBarriers.at(i1)/p.get_T_dust()),
                   s.quantMobilities.at(i1)))
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density())
       / p.get_dust2gas_num() * r.abc[0] * calc_cross_surf_barrier_prob(p, r);
  return NV_Ith_S(y, i0)
       * NV_Ith_S(y, i1)
       * tmp;
}


std::vector<TYPES::DTP_FLOAT> drdy_surface_AB(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  int i0 = r.idxReactants[0], i1 = r.idxReactants[1];
  TYPES::DTP_FLOAT tmp =
         (std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/p.get_T_dust()),
                   s.quantMobilities.at(i0))
        + std::max(s.vibFreqs.at(i1) * exp(-s.diffBarriers.at(i1)/p.get_T_dust()),
                   s.quantMobilities.at(i1)))
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density())
       / p.get_dust2gas_num() * r.abc[0] * calc_cross_surf_barrier_prob(p, r);
  drdy.push_back(tmp * NV_Ith_S(y, i1));
  drdy.push_back(tmp * NV_Ith_S(y, i0));
  return drdy;
}


TYPES::DTP_FLOAT rate_surf2mant(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  int i0 = r.idxReactants[0];
  if (m.t_calc != t) {
    m.t_calc = t;
    m.k_ads_tot = 0.0;
    m.k_eva_tot = 0.0;
    m.surf_tot = 0.0;
    m.mant_tot = 0.0;
    for (auto const& a: m.ads_reactions) {
      m.k_ads_tot += rate_adsorption(t, y, a, p, s, m);
    }
    for (auto const& a: m.eva_reactions) {
      m.k_eva_tot += rate_desorption(t, y, a, p, s, m);
    }
    for (auto const& i: s.surfaceSpecies) {
      m.surf_tot += NV_Ith_S(y, i);
    }
    for (auto const& i: s.mantleSpecies) {
      m.mant_tot += NV_Ith_S(y, i);
    }
    std::cout << t << " " << m.k_ads_tot << " " << m.mant_tot << std::endl;
  }
  return NV_Ith_S(y, i0)
       * m.k_ads_tot / p.get_dust2gas_num()
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density());
}


std::vector<TYPES::DTP_FLOAT> drdy_surf2mant(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(
         m.k_ads_tot / p.get_dust2gas_num()
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density()));
  return drdy;
}


TYPES::DTP_FLOAT rate_mant2surf(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  int i0 = r.idxReactants[0];
  if (m.t_calc != t) {
    m.t_calc = t;
    m.k_ads_tot = 0.0;
    m.k_eva_tot = 0.0;
    m.surf_tot = 0.0;
    m.mant_tot = 0.0;
    for (auto const& a: m.ads_reactions) {
      m.k_ads_tot += rate_adsorption(t, y, a, p, s, m);
    }
    for (auto const& a: m.eva_reactions) {
      m.k_eva_tot += rate_desorption(t, y, a, p, s, m);
    }
    for (auto const& i: s.surfaceSpecies) {
      m.surf_tot += NV_Ith_S(y, i);
    }
    for (auto const& i: s.mantleSpecies) {
      m.mant_tot += NV_Ith_S(y, i);
    }
  }
  return NV_Ith_S(y, i0)
       * m.k_eva_tot / p.get_dust2gas_num()
       / (std::max(m.surf_tot, m.mant_tot) + 1e-40);
}


std::vector<TYPES::DTP_FLOAT> drdy_mant2surf(
    const TYPES::DTP_FLOAT& t,
    const TYPES::DTP_Y& y,
    const TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  std::vector<TYPES::DTP_FLOAT> drdy;
  drdy.push_back(
         m.k_eva_tot / p.get_dust2gas_num()
       / (std::max(m.surf_tot, m.mant_tot) + 1e-40));
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
  (*user_data.rate_calculators)[5]  = rate_ion_neutral;
  (*user_data.drdy_calculators)[5]  = drdy_ion_neutral;
  (*user_data.rate_calculators)[61] = rate_adsorption;
  (*user_data.drdy_calculators)[61] = drdy_adsorption;
  (*user_data.rate_calculators)[62] = rate_desorption;
  (*user_data.drdy_calculators)[62] = drdy_desorption;
  (*user_data.rate_calculators)[63] = rate_surface_AA;
  (*user_data.drdy_calculators)[63] = drdy_surface_AA;
  (*user_data.rate_calculators)[64] = rate_surface_AB;
  (*user_data.drdy_calculators)[64] = drdy_surface_AB;
  (*user_data.rate_calculators)[65] = rate_mant2surf;
  (*user_data.drdy_calculators)[65] = drdy_mant2surf;
  (*user_data.rate_calculators)[66] = rate_surf2mant;
  (*user_data.drdy_calculators)[66] = drdy_surf2mant;
  //(*user_data.rate_calculators)[3]  = rate_photoionization;
  //(*user_data.drdy_calculators)[3]  = drdy_photoionization;
  //(*user_data.rate_calculators)[73] = rate_photoionization;
  //(*user_data.drdy_calculators)[73] = drdy_photoionization;
  //(*user_data.rate_calculators)[53] = rate_ion_neutral;
  //(*user_data.drdy_calculators)[53] = drdy_ion_neutral;
  return 0;
}

}

#endif //CALC_RATE_H
