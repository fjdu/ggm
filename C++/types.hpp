#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include "constants.hpp"

namespace TYPES {

typedef double DTP_FLOAT;
typedef double* DTP_Y;


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

  PhyParams& set_n_gas(DTP_FLOAT _v);
  PhyParams& set_T_gas(DTP_FLOAT _v);
  PhyParams& set_T_dust(DTP_FLOAT _v);
  PhyParams& set_Av(DTP_FLOAT _v);
  PhyParams& set_G0_UV(DTP_FLOAT _v);
  PhyParams& set_chi_Xray(DTP_FLOAT _v);
  PhyParams& set_chi_cosmicray(DTP_FLOAT _v);
  PhyParams& set_dust2gas_num(DTP_FLOAT _v);
  PhyParams& set_dust2gas_mass(DTP_FLOAT _v);
  PhyParams& set_dust_material_density(DTP_FLOAT _v);
  PhyParams& set_dust_site_density(DTP_FLOAT _v);
  PhyParams& set_dust_radius(DTP_FLOAT _v);
  PhyParams& set_dust_crosssec(DTP_FLOAT _v);
  PhyParams& set_dust_albedo(DTP_FLOAT _v);
  PhyParams& set_mean_mol_weight(DTP_FLOAT _v);
  PhyParams& set_chemdesorption_factor(DTP_FLOAT _v);

  void prep_params();

  int from_file(std::string fname);

private:
  DTP_FLOAT n_gas, T_gas, T_dust, Av, G0_UV, chi_Xray, chi_cosmicray,
    dust2gas_num, dust2gas_mass, dust_material_density,
    dust_site_density, dust_radius, dust_crosssec, dust_albedo,
    mean_mol_weight, chemdesorption_factor;
};


extern std::map<std::string, void (*)(PhyParams&, DTP_FLOAT)> phySetterDict;
extern std::map<std::string, DTP_FLOAT (*)(PhyParams&)> phyGetterDict;


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
    double drdy[2];
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
    Species() {
    }
    ~Species() {
    }
};


class OtherData {
  public:
    DTP_FLOAT t_calc, k_eva_tot, k_ads_tot, mant_tot, surf_tot;
    Reactions ads_reactions, eva_reactions;
};


typedef DTP_FLOAT (*RateCalculator)(
    const DTP_FLOAT&, double *, Reaction&,
    const PhyParams&, const Species&, TYPES::OtherData&);
typedef std::vector<DTP_FLOAT> (*dRdyCalculator)(
    const DTP_FLOAT&, double *, Reaction&,
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
