#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include "constants.hpp"
#include "types.hpp"

namespace TYPES {

PhyParams& PhyParams::set_n_gas(DTP_FLOAT _v) {n_gas = _v; return *this;}
PhyParams& PhyParams::set_T_gas(DTP_FLOAT _v) {T_gas = _v; return *this;}
PhyParams& PhyParams::set_T_dust(DTP_FLOAT _v) {T_dust = _v; return *this;}
PhyParams& PhyParams::set_Av(DTP_FLOAT _v) {Av = _v; return *this;}
PhyParams& PhyParams::set_G0_UV(DTP_FLOAT _v) {G0_UV = _v; return *this;}
PhyParams& PhyParams::set_chi_Xray(DTP_FLOAT _v) {chi_Xray = _v; return *this;}
PhyParams& PhyParams::set_chi_cosmicray(DTP_FLOAT _v) {chi_cosmicray = _v; return *this;}
PhyParams& PhyParams::set_dust2gas_num(DTP_FLOAT _v) {dust2gas_num = _v; return *this;}
PhyParams& PhyParams::set_dust2gas_mass(DTP_FLOAT _v) {dust2gas_mass = _v; return *this;}
PhyParams& PhyParams::set_dust_material_density(DTP_FLOAT _v) {dust_material_density = _v; return *this;}
PhyParams& PhyParams::set_dust_site_density(DTP_FLOAT _v) {dust_site_density = _v; return *this;}
PhyParams& PhyParams::set_dust_radius(DTP_FLOAT _v) {dust_radius = _v; return *this;}
PhyParams& PhyParams::set_dust_crosssec(DTP_FLOAT _v) {dust_crosssec = _v; return *this;}
PhyParams& PhyParams::set_dust_albedo(DTP_FLOAT _v) {dust_albedo = _v; return *this;}
PhyParams& PhyParams::set_mean_mol_weight(DTP_FLOAT _v) {mean_mol_weight = _v; return *this;}
PhyParams& PhyParams::set_chemdesorption_factor(DTP_FLOAT _v) {chemdesorption_factor = _v; return *this;}
PhyParams& PhyParams::set_t_max_year(DTP_FLOAT _v) {t_max_year = _v; return *this;}


DTP_FLOAT PhyParams::get_T_gas(DTP_FLOAT t) const {
  return T_gas;
  double t_scale = CONST::phy_SecondsPerYear * 1e5;
  double ts = t / t_scale;
  double Tbase = 15.0, Tscale=70.0;
  double T = Tbase + Tscale * ts * ts * exp(-ts);
  return T;
}


DTP_FLOAT PhyParams::get_T_dust(DTP_FLOAT t) const {
  return T_dust;
  double t_scale = CONST::phy_SecondsPerYear * 1e5;
  double ts = t / t_scale;
  double Tbase = 10.0, Tscale=70.0;
  double T = Tbase + Tscale * ts * ts * exp(-ts);
  return T;
}


DTP_FLOAT PhyParams::get_n_gas(DTP_FLOAT t) const {
  return n_gas;
  double t_scale = CONST::phy_SecondsPerYear * 1e5;
  double ts = t / t_scale;
  double nbase = 1e6;
  double n = nbase * (1.0 + ts);
  return n;
}


inline DTP_FLOAT get_n_gas(PhyParams& p) {return p.get_n_gas();}
inline DTP_FLOAT get_T_gas(PhyParams& p) {return p.get_T_gas();}
inline DTP_FLOAT get_T_dust(PhyParams& p) {return p.get_T_dust();}
inline DTP_FLOAT get_Av(PhyParams& p) {return p.get_Av();}
inline DTP_FLOAT get_G0_UV(PhyParams& p) {return p.get_G0_UV();}
inline DTP_FLOAT get_chi_Xray(PhyParams& p) {return p.get_chi_Xray();}
inline DTP_FLOAT get_chi_cosmicray(PhyParams& p) {return p.get_chi_cosmicray();}
inline DTP_FLOAT get_dust2gas_num(PhyParams& p) {return p.get_dust2gas_num();}
inline DTP_FLOAT get_dust2gas_mass(PhyParams& p) {return p.get_dust2gas_mass();}
inline DTP_FLOAT get_dust_material_density(PhyParams& p) {return p.get_dust_material_density();}
inline DTP_FLOAT get_dust_site_density(PhyParams& p) {return p.get_dust_site_density();}
inline DTP_FLOAT get_dust_radius(PhyParams& p) {return p.get_dust_radius();}
inline DTP_FLOAT get_dust_crosssec(PhyParams& p) {return p.get_dust_crosssec();}
inline DTP_FLOAT get_dust_albedo(PhyParams& p) {return p.get_dust_albedo();}
inline DTP_FLOAT get_mean_mol_weight(PhyParams& p) {return p.get_mean_mol_weight();}
inline DTP_FLOAT get_chemdesorption_factor(PhyParams& p) {return p.get_chemdesorption_factor();}
inline DTP_FLOAT get_t_max_year(PhyParams& p) {return p.get_t_max_year();}


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
void set_t_max_year(PhyParams& p, DTP_FLOAT _v) {p.set_t_max_year(_v);}


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
  {"chemdesorption_factor", set_chemdesorption_factor},
  {"t_max_year", set_t_max_year}
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
  {"chemdesorption_factor", get_chemdesorption_factor},
  {"t_max_year", get_t_max_year}
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
    std::cout << "Error in PhyParams::from_file: " << fname << std::endl;
  }
  return 0;
}

}
