#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include "sundials/sundials_types.h"
#include "sundials/sundials_config.h"
#include "sunmatrix/sunmatrix_sparse.h"
#include "nvector/nvector_serial.h"
#include "constants.hpp"

namespace TYPES {

typedef double DTP_FLOAT;
typedef N_Vector DTP_Y;
typedef SUNMatrix DTP_Jac;


class DPT_Config {
  std::string f_reactions, f_initial_abundances;
};


class PhyParams {
public:
  inline DTP_FLOAT get_n_gas() const {return n_gas;}
  inline DTP_FLOAT get_T_gas() const {return T_gas;}
  inline DTP_FLOAT get_T_dust() const {return T_dust;}
  inline DTP_FLOAT get_Av() const {return Av;}
  inline DTP_FLOAT get_G0_UV() const {return G0_UV;}
  inline DTP_FLOAT get_chi_Xray() const {return chi_Xray;}
  inline DTP_FLOAT get_chi_cosmicray() const {return chi_cosmicray;}
  inline DTP_FLOAT get_dust2gas_num() const {return dust2gas_num;}
  inline DTP_FLOAT get_dust2gas_mass() const {return dust2gas_mass;}
  inline DTP_FLOAT get_dust_material_density() const {return dust_material_density;}
  inline DTP_FLOAT get_dust_site_density() const {return dust_site_density;}
  inline DTP_FLOAT get_dust_radius() const {return dust_radius;}
  inline DTP_FLOAT get_dust_crosssec() const {return dust_crosssec;}
  inline DTP_FLOAT get_dust_albedo() const {return dust_albedo;}
  inline DTP_FLOAT get_mean_mol_weight() const {return mean_mol_weight;}
  inline DTP_FLOAT get_chemdesorption_factor() const {return chemdesorption_factor;}

  inline DTP_FLOAT get_n_gas(DTP_FLOAT t) const;
  inline DTP_FLOAT get_T_gas(DTP_FLOAT t) const;
  inline DTP_FLOAT get_T_dust(DTP_FLOAT t) const;
  inline DTP_FLOAT get_Av(DTP_FLOAT t) const;
  inline DTP_FLOAT get_G0_UV(DTP_FLOAT t) const;
  inline DTP_FLOAT get_chi_Xray(DTP_FLOAT t) const;
  inline DTP_FLOAT get_chi_cosmicray(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust2gas_num(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust2gas_mass(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_material_density(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_site_density(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_radius(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_crosssec(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_albedo(DTP_FLOAT t) const;
  inline DTP_FLOAT get_mean_mol_weight(DTP_FLOAT t) const;
  inline DTP_FLOAT get_chemdesorption_factor(DTP_FLOAT t) const;

  PhyParams& set_n_gas(DTP_FLOAT _v) {n_gas = _v; return *this;}
  PhyParams& set_T_gas(DTP_FLOAT _v) {T_gas = _v; return *this;}
  PhyParams& set_T_dust(DTP_FLOAT _v) {T_dust = _v; return *this;}
  PhyParams& set_Av(DTP_FLOAT _v) {Av = _v; return *this;}
  PhyParams& set_G0_UV(DTP_FLOAT _v) {G0_UV = _v; return *this;}
  PhyParams& set_chi_Xray(DTP_FLOAT _v) {chi_Xray = _v; return *this;}
  PhyParams& set_chi_cosmicray(DTP_FLOAT _v) {chi_cosmicray = _v; return *this;}
  PhyParams& set_dust2gas_num(DTP_FLOAT _v) {dust2gas_num = _v; return *this;}
  PhyParams& set_dust2gas_mass(DTP_FLOAT _v) {dust2gas_mass = _v; return *this;}
  PhyParams& set_dust_material_density(DTP_FLOAT _v) {dust_material_density = _v; return *this;}
  PhyParams& set_dust_site_density(DTP_FLOAT _v) {dust_site_density = _v; return *this;}
  PhyParams& set_dust_radius(DTP_FLOAT _v) {dust_radius = _v; return *this;}
  PhyParams& set_dust_crosssec(DTP_FLOAT _v) {dust_crosssec = _v; return *this;}
  PhyParams& set_dust_albedo(DTP_FLOAT _v) {dust_albedo = _v; return *this;}
  PhyParams& set_mean_mol_weight(DTP_FLOAT _v) {mean_mol_weight = _v; return *this;}
  PhyParams& set_chemdesorption_factor(DTP_FLOAT _v) {chemdesorption_factor = _v; return *this;}

  void prep_params();

  int from_file(std::string fname);

private:
  DTP_FLOAT n_gas, T_gas, T_dust, Av, G0_UV, chi_Xray, chi_cosmicray,
    dust2gas_num, dust2gas_mass, dust_material_density,
    dust_site_density, dust_radius, dust_crosssec, dust_albedo,
    mean_mol_weight, chemdesorption_factor;
};


void set_n_gas(PhyParams& p, DTP_FLOAT _v) {p.set_n_gas(_v);}
void set_T_gas(PhyParams& p, DTP_FLOAT _v) {p.set_T_gas(_v);}
void set_T_dust(PhyParams& p, DTP_FLOAT _v) {p.set_T_dust(_v);}
void set_Av(PhyParams& p, DTP_FLOAT _v) {p.set_Av(_v);}
void set_G0_UV(PhyParams& p, DTP_FLOAT _v) {p.set_G0_UV(_v);}
void set_chi_Xray(PhyParams& p, DTP_FLOAT _v) {p.set_chi_Xray(_v);}
void set_chi_cosmicray(PhyParams& p, DTP_FLOAT _v) {p.set_chi_cosmicray(_v);}
void set_dust2gas_num(PhyParams& p, DTP_FLOAT _v) {p.set_dust2gas_num(_v);}
void set_dust2gas_mass(PhyParams& p, DTP_FLOAT _v) {p.set_dust2gas_mass(_v);}
void set_dust_material_density(PhyParams& p, DTP_FLOAT _v) {p.set_dust_material_density(_v);}
void set_dust_site_density(PhyParams& p, DTP_FLOAT _v) {p.set_dust_site_density(_v);}
void set_dust_radius(PhyParams& p, DTP_FLOAT _v) {p.set_dust_radius(_v);}
void set_dust_crosssec(PhyParams& p, DTP_FLOAT _v) {p.set_dust_crosssec(_v);}
void set_dust_albedo(PhyParams& p, DTP_FLOAT _v) {p.set_dust_albedo(_v);}
void set_mean_mol_weight(PhyParams& p, DTP_FLOAT _v) {p.set_mean_mol_weight(_v);}
void set_chemdesorption_factor(PhyParams& p, DTP_FLOAT _v) {p.set_chemdesorption_factor(_v);}


DTP_FLOAT get_n_gas(PhyParams& p) {return p.get_n_gas();}
DTP_FLOAT get_T_gas(PhyParams& p) {return p.get_T_gas();}
DTP_FLOAT get_T_dust(PhyParams& p) {return p.get_T_dust();}
DTP_FLOAT get_Av(PhyParams& p) {return p.get_Av();}
DTP_FLOAT get_G0_UV(PhyParams& p) {return p.get_G0_UV();}
DTP_FLOAT get_chi_Xray(PhyParams& p) {return p.get_chi_Xray();}
DTP_FLOAT get_chi_cosmicray(PhyParams& p) {return p.get_chi_cosmicray();}
DTP_FLOAT get_dust2gas_num(PhyParams& p) {return p.get_dust2gas_num();}
DTP_FLOAT get_dust2gas_mass(PhyParams& p) {return p.get_dust2gas_mass();}
DTP_FLOAT get_dust_material_density(PhyParams& p) {return p.get_dust_material_density();}
DTP_FLOAT get_dust_site_density(PhyParams& p) {return p.get_dust_site_density();}
DTP_FLOAT get_dust_radius(PhyParams& p) {return p.get_dust_radius();}
DTP_FLOAT get_dust_crosssec(PhyParams& p) {return p.get_dust_crosssec();}
DTP_FLOAT get_dust_albedo(PhyParams& p) {return p.get_dust_albedo();}
DTP_FLOAT get_mean_mol_weight(PhyParams& p) {return p.get_mean_mol_weight();}
DTP_FLOAT get_chemdesorption_factor(PhyParams& p) {return p.get_chemdesorption_factor();}


std::map<std::string, void (*)(PhyParams&, DTP_FLOAT)> phySetterDict = {
  {"n_gas", set_n_gas},
  {"T_gas", set_T_gas},
  {"T_dust", set_T_dust},
  {"Av", set_Av},
  {"G0_UV", set_G0_UV},
  {"chi_Xray", set_chi_Xray},
  {"chi_cosmicray", set_chi_cosmicray},
  {"dust2gas_num", set_dust2gas_num},
  {"dust2gas_mass", set_dust2gas_mass},
  {"dust_material_density", set_dust_material_density},
  {"dust_site_density", set_dust_site_density},
  {"dust_radius", set_dust_radius},
  {"dust_crosssec", set_dust_crosssec},
  {"dust_albedo", set_dust_albedo},
  {"mean_mol_weight", set_mean_mol_weight},
  {"chemdesorption_factor", set_chemdesorption_factor}
};


std::map<std::string, DTP_FLOAT (*)(PhyParams&)> phyGetterDict = {
  {"n_gas", get_n_gas},
  {"T_gas", get_T_gas},
  {"T_dust", get_T_dust},
  {"Av", get_Av},
  {"G0_UV", get_G0_UV},
  {"chi_Xray", get_chi_Xray},
  {"chi_cosmicray", get_chi_cosmicray},
  {"dust2gas_num", get_dust2gas_num},
  {"dust2gas_mass", get_dust2gas_mass},
  {"dust_material_density", get_dust_material_density},
  {"dust_site_density", get_dust_site_density},
  {"dust_radius", get_dust_radius},
  {"dust_crosssec", get_dust_crosssec},
  {"dust_albedo", get_dust_albedo},
  {"mean_mol_weight", get_mean_mol_weight},
  {"chemdesorption_factor", get_chemdesorption_factor}
};


void PhyParams::prep_params() {
  dust2gas_num = dust2gas_mass * (CONST::phy_mProton_CGS * mean_mol_weight)
              / (4.0*CONST::PI/3.0 * dust_material_density *
                 dust_radius * dust_radius * dust_radius);
  dust_crosssec = CONST::PI * dust_radius * dust_radius;
}


int PhyParams::from_file(std::string fname)
{
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s*(=)\s*(\S+))");
  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {
      if (std::regex_match(line, comment)) {
        continue;
      }
      if (std::regex_match(line, emptyline)) {
        continue;
      }
      std::smatch matches;
      if (std::regex_search(line, matches, entry)) {
        std::string var_name = matches.str(1);
        if (phySetterDict.find(var_name) != phySetterDict.end()) {
          DTP_FLOAT var_val = std::stod(matches.str(3));
          phySetterDict[var_name](*this, var_val);
        } else {
          std::cout << "Invalid line: " << line << std::endl;
        }
      }
    }
  } else {
    std::cout << "Error reading file: " << fname << std::endl;
  }
  return 0;
}


class ReactionFormat {
  public:
    int nReactants, nProducts, nABC, lenSpeciesName, rowlen_min;
    ReactionFormat(int nReactants=3, int nProducts=4, int nABC=3,
        int lenSpeciesName=12, int lenABC=9, int lenT=6,
        int rowlen_min=126)
    {
    }
};


class Reaction {
  public:
    int itype;
    std::vector<int> idxReactants, idxProducts;
    std::vector<DTP_FLOAT> abc;
    std::vector<DTP_FLOAT> Trange;
};


typedef std::vector<Reaction> Reactions;
typedef std::map<int, int> ReactionTypes;
typedef std::map<std::string, DTP_FLOAT> Elements;


class Species {
  public:
    std::map<std::string, int> name2idx;
    std::map<int, std::string> idx2name;
    std::map<int, std::map<std::string, int> > elementsSpecies;
    std::map<int, DTP_FLOAT> massSpecies, enthalpies,
                             vibFreqs, diffBarriers, quantMobilities;
    std::set<int> gasSpecies, surfaceSpecies, mantleSpecies;
    std::vector<DTP_FLOAT> abundances;
    DTP_Y y;
};


class OtherData {
  public:
    DTP_FLOAT t_calc, k_eva_tot, k_ads_tot, mant_tot, surf_tot;
    Reactions ads_reactions, eva_reactions;
};


typedef DTP_FLOAT (*RateCalculator)(
    const DTP_FLOAT&, const DTP_Y&, const Reaction&,
    const PhyParams&, const Species&, TYPES::OtherData&);
typedef std::vector<DTP_FLOAT> (*dRdyCalculator)(
    const DTP_FLOAT&, const DTP_Y&, const Reaction&,
    const PhyParams&, const Species&, TYPES::OtherData&);


typedef std::map<int, RateCalculator> RateCalculators;
typedef std::map<int, dRdyCalculator> dRdyCalculators;


class User_data {
  public:
    PhyParams* physical_params;
    Reactions* reactions;
    ReactionTypes* reaction_types;
    Species* species;
    RateCalculators* rate_calculators;
    dRdyCalculators* drdy_calculators;
    OtherData* other_data;

    User_data() {
      physical_params = new PhyParams();
      reactions = new Reactions();
      reaction_types = new ReactionTypes();
      species = new Species();
      rate_calculators = new RateCalculators();
      drdy_calculators = new dRdyCalculators();
      other_data = new OtherData();
    }
    ~User_data() {
      // delete stuff
      delete physical_params;
      delete reactions;
      delete reaction_types;
      delete species;
      delete rate_calculators;
      delete drdy_calculators;
      delete other_data;
    }
  private:
};


}
#endif //TYPES_H
