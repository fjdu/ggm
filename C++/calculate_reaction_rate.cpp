#include "math.h"
#include <algorithm>
#include "types.hpp"
#include "constants.hpp"

#define NV_Ith_S(x, i) (x[i])

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
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  r.drdy[0] = r.abc[0] // Act as sticking efficiency
       * p.get_n_gas(t)
       * p.get_dust2gas_num()
       * p.get_dust_crosssec()
       * thermal_velocity_CGS(p.get_T_gas(t),
                              s.massSpecies.at(r.idxReactants[0]));
  return NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
}


TYPES::DTP_FLOAT rate_desorption(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  r.drdy[0] = s.vibFreqs.at(r.idxReactants[0])
    * r.abc[0]
    * (exp(-r.abc[2]/p.get_T_dust(t))
     + p.get_chi_cosmicray()
     * CONST::phy_cosmicray_desorption_factor
     * exp(-r.abc[2] / CONST::phy_cosmicray_desorption_T));
  return r.drdy[0] * NV_Ith_S(y, r.idxReactants[0]);
}


TYPES::DTP_FLOAT rate_ion_neutral(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  double Tgas = p.get_T_gas(t), tmp;
  if ((Tgas < r.Trange[0]) || (Tgas > r.Trange[1])) {
    tmp = 0.0;
  } else {
    tmp = p.get_n_gas(t) * r.abc[0]
        * pow(p.get_T_gas(t)/3e2, r.abc[1])
        * exp(-r.abc[2] / p.get_T_gas(t));
  }
  r.drdy[0] = tmp * NV_Ith_S(y, r.idxReactants[1]);
  r.drdy[1] = tmp * NV_Ith_S(y, r.idxReactants[0]);
  return NV_Ith_S(y, r.idxReactants[0])
       * NV_Ith_S(y, r.idxReactants[1])
       * tmp;
}


TYPES::DTP_FLOAT rate_cosmicray_ionization(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  r.drdy[0] = r.abc[0] * p.get_chi_cosmicray();
  return NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
}


TYPES::DTP_FLOAT rate_cosmicray_induced_ionization(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  r.drdy[0] = r.abc[0] * p.get_chi_cosmicray()
       * pow(p.get_T_gas(t)/3e2, r.abc[1])
       * r.abc[2] / (1.0 - p.get_dust_albedo());
  return NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
}


TYPES::DTP_FLOAT rate_photoionization(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  r.drdy[0] = r.abc[0] * p.get_G0_UV() * exp(-r.abc[2] * p.get_Av());
  return NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
}


inline TYPES::DTP_FLOAT calc_cross_surf_barrier_prob(
    const double Tdust, TYPES::Reaction& r)
{
  TYPES::DTP_FLOAT tmp = 1.0;
  if ((r.abc[2] > 0.0) && (r.Trange[0] > 0.0)) {
    tmp =
      exp(std::max(-r.abc[2]/Tdust,
                   -2.0 * r.abc[1] * CONST::phy_DiffBarrierWidth_CGS
                   / CONST::phy_hbarPlanck_CGS
                   * sqrt(2.0 * r.Trange[0] * CONST::phy_mProton_CGS
                        * CONST::phy_kBoltzmann_CGS * r.abc[2])));
  }
  return tmp;
}


TYPES::DTP_FLOAT rate_surface_AA(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  double Tdust = p.get_T_dust(t);
  int i0 = r.idxReactants[0];
  TYPES::DTP_FLOAT tmp =
         std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/Tdust),
                  s.quantMobilities.at(i0))
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density())
       / p.get_dust2gas_num() * r.abc[0] * calc_cross_surf_barrier_prob(Tdust, r);
  r.drdy[0] = tmp * NV_Ith_S(y, i0);
  r.drdy[1] = tmp * NV_Ith_S(y, i0);
  return NV_Ith_S(y, i0) * NV_Ith_S(y, i0) * tmp;
}


TYPES::DTP_FLOAT rate_surface_AB(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  double Tdust = p.get_T_dust(t);
  int i0 = r.idxReactants[0], i1 = r.idxReactants[1];
  TYPES::DTP_FLOAT tmp =
      (std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/Tdust),
                s.quantMobilities.at(i0))
     + std::max(s.vibFreqs.at(i1) * exp(-s.diffBarriers.at(i1)/Tdust),
                s.quantMobilities.at(i1)))
    / (4.0*p.get_dust_crosssec()*p.get_dust_site_density())
    / p.get_dust2gas_num() * r.abc[0] * calc_cross_surf_barrier_prob(Tdust, r);
  r.drdy[0] = tmp * NV_Ith_S(y, i1);
  r.drdy[1] = tmp * NV_Ith_S(y, i0);
  return NV_Ith_S(y, i0) * NV_Ith_S(y, i1) * tmp;
}


void update_surfmant(
    const TYPES::DTP_FLOAT& t,
    double *y,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m) {
  m.t_calc = t;
  m.k_ads_tot = 0.0;
  m.k_eva_tot = 0.0;
  m.surf_tot = 0.0;
  m.mant_tot = 0.0;
  for (auto & a: m.ads_reactions) {
    m.k_ads_tot += rate_adsorption(t, y, a, p, s, m);
  }
  for (auto & a: m.eva_reactions) {
    m.k_eva_tot += rate_desorption(t, y, a, p, s, m);
  }
  for (auto const& i: s.surfaceSpecies) {
    m.surf_tot += NV_Ith_S(y, i);
  }
  for (auto const& i: s.mantleSpecies) {
    m.mant_tot += NV_Ith_S(y, i);
  }
}


TYPES::DTP_FLOAT rate_surf2mant(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  int i0 = r.idxReactants[0];
  if (m.t_calc != t) {
    update_surfmant(t, y, p, s, m);
  }
  r.drdy[0] = m.k_ads_tot / p.get_dust2gas_num()
       / (4.0*p.get_dust_crosssec()*p.get_dust_site_density());
  return NV_Ith_S(y, i0) * r.drdy[0];
}


TYPES::DTP_FLOAT rate_mant2surf(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  int i0 = r.idxReactants[0];
  if (m.t_calc != t) {
    update_surfmant(t, y, p, s, m);
  }
  r.drdy[0] = m.k_eva_tot
       / (std::max(m.surf_tot, m.mant_tot) + 1e-60);
  return NV_Ith_S(y, i0) * r.drdy[0];
}


TYPES::DTP_FLOAT rate_iongrain(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  double TemperatureReduced = CONST::phy_kBoltzmann_SI * p.get_T_gas(t)
                            / (pow(CONST::phy_elementaryCharge_SI,2)
                             * CONST::phy_CoulombConst_SI
                             / (p.get_dust_radius()*1e-2));
  double JNegaPosi = (1.0 + 1.0/TemperatureReduced)
                   * (1.0 + sqrt(2.0/(2.0+TemperatureReduced)));
  double JChargeNeut = 1.0 + sqrt(CONST::PI/2.0/TemperatureReduced);

  int i0 = r.idxReactants[0], i1 = r.idxReactants[1], i2;
  int charge = s.elementsSpecies.at(i0).at("+") * s.elementsSpecies.at(i1).at("-") +
               s.elementsSpecies.at(i0).at("-") * s.elementsSpecies.at(i1).at("+");
  if (s.elementsSpecies.at(i0).at("Grain") == 0) {
    i2 = i0;
  } else if (s.elementsSpecies.at(i1).at("Grain") == 0) {
    i2 = i1;
  } else {
    std::cout << "Error in rate_iongrain." << std::endl;
  }
  int tmp;
  if (charge < 0) {
    tmp = p.get_dust_crosssec() * JNegaPosi
        * thermal_velocity_CGS(p.get_T_gas(t), s.massSpecies.at(i2));
  } else {
    tmp = p.get_dust_crosssec() * JChargeNeut
        * thermal_velocity_CGS(p.get_T_gas(t), s.massSpecies.at(i2));
  }
  r.drdy[0] = tmp * NV_Ith_S(y, i1);
  r.drdy[1] = tmp * NV_Ith_S(y, i0);
  return NV_Ith_S(y, i0) * NV_Ith_S(y, i1) * tmp;
}


TYPES::DTP_FLOAT rate_photodesorption(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  double yield = r.abc[0] + r.abc[1] * p.get_T_gas(t);
  double uv_flux = p.get_G0_UV() * exp(-CONST::phy_UVext2Av * p.get_Av())
                 * CONST::phy_Habing_photon_flux_CGS;
  double rate = uv_flux * yield * p.get_dust2gas_num() * p.get_dust_crosssec();
  double tmp1 = p.get_dust2gas_num()
              * (4.0 * p.get_dust_site_density() * p.get_dust_crosssec());
  double tmp = NV_Ith_S(y, r.idxReactants[0]) / tmp1;
  if (tmp < 1e-6) {
    r.drdy[0] = rate / tmp1;
    return rate * tmp;
  } else {
    double tmp2 = exp(-tmp);
    r.drdy[0] = rate * tmp2 / tmp1;
    return rate * (1.0 - tmp2);
  }
}


TYPES::DTP_FLOAT rate_dummy(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::OtherData& m)
{
  r.drdy[0] = 0.0;
  r.drdy[1] = 0.0;
  return 0.0;
}


int assignReactionHandlers(TYPES::User_data& user_data) {
  (*user_data.rate_calculators)[1]  = rate_cosmicray_ionization;
  (*user_data.rate_calculators)[71] = rate_cosmicray_ionization;
  (*user_data.rate_calculators)[2]  = rate_cosmicray_induced_ionization;
  (*user_data.rate_calculators)[72] = rate_cosmicray_induced_ionization;
  (*user_data.rate_calculators)[5]  = rate_ion_neutral;
  (*user_data.rate_calculators)[61] = rate_adsorption;
  (*user_data.rate_calculators)[62] = rate_desorption;
  (*user_data.rate_calculators)[63] = rate_surface_AA;
  (*user_data.rate_calculators)[64] = rate_surface_AB;
  (*user_data.rate_calculators)[65] = rate_mant2surf;
  (*user_data.rate_calculators)[66] = rate_surf2mant;
  (*user_data.rate_calculators)[3]  = rate_photoionization;
  (*user_data.rate_calculators)[20] = rate_cosmicray_induced_ionization;
  (*user_data.rate_calculators)[21] = rate_iongrain;
  //(*user_data.rate_calculators)[53] = rate_ion_neutral;
  (*user_data.rate_calculators)[75] = rate_photodesorption;
  //
  (*user_data.rate_calculators)[53] = rate_dummy;
  //(*user_data.rate_calculators)[65] = rate_dummy;
  //(*user_data.rate_calculators)[66] = rate_dummy;
  (*user_data.rate_calculators)[13] = rate_dummy;
  (*user_data.rate_calculators)[67] = rate_dummy;
  return 0;
}

}
