#include "logistics.cpp"
#include "types.hpp"
#include "rate_equation.cpp"
#include "constants.hpp"
#include "calculate_reaction_rate.cpp"
#include <iostream>
#include <string>
#include <algorithm>

//DYLD_LIBRARY_PATH=/Users/fjdu/Dropbox/codes_from_others/sundials-4.1.0/builddir/lib/ ./re

int main(int argc, char **argv)
{
  for (int i=0; i<argc; ++i) {
    std::cout << std::string(argv[i]) << " ";
  }
  std::cout << std::endl;

  TYPES::User_data user_data;

  user_data.physical_params->from_file("phy_params.dat");
  user_data.physical_params->prep_params();
  for (auto const& p: TYPES::phySetterDict) {
    std::cout << p.first << " = " <<
      TYPES::phyGetterDict[p.first](*user_data.physical_params) << " " << std::endl;
  }

  LOGIS::load_reactions("/Users/fjdu/_c/DeuterationGasGrain/output_1/rate06_dipole_reduced_20120513_50K_com_red_com.dat", &user_data);
  std::cout << "Number of reactions: "
            << user_data.reactions->size() << std::endl;
  std::cout << "Number of species: "
            << user_data.species->name2idx.size() << std::endl;
  std::cout << "Number of reaction types: "
            << user_data.reaction_types->size() << std::endl;
  //for (auto const& i: *(user_data.reaction_types)) {
  //  std::cout << i << " ";
  //}
  //std::cout << std::endl;

  LOGIS::assignElementsToSpecies(*(user_data.species), CONST::element_masses);
  LOGIS::calculateSpeciesMasses(*(user_data.species), CONST::element_masses);
  LOGIS::calculateSpeciesVibFreqs(*(user_data.species), *(user_data.reactions));
  //for (auto const& s: user_data.species->name2idx) {
  //  std::cout << s.first << " ";
  //  auto eleDict = user_data.species->elementsSpecies[s.second];
  //  for (auto const& e: CONST::element_masses) {
  //    if (eleDict[e.first] != 0) {
  //      std::cout << e.first << "(" << eleDict[e.first] << ") ";
  //    }
  //  }
  //  std::cout << " " << user_data.species->massSpecies[s.second] << std::endl;
  //}

  LOGIS::classifySpeciesByPhase(*(user_data.species));
  std::cout << "Number of gas species: "
            << user_data.species->gasSpecies.size() << std::endl;
  std::cout << "Number of surface species: "
            << user_data.species->surfaceSpecies.size() << std::endl;
  std::cout << "Number of mantle species: "
            << user_data.species->mantleSpecies.size() << std::endl;

  LOGIS::loadInitialAbundances(*user_data.species, "initial_abundances.dat");
  std::cout << "Number of species with initial abundances: "
            << std::count_if(user_data.species->abundances.begin(),
                             user_data.species->abundances.end(),
                             [](TYPES::DTP_FLOAT v){return v>1e-40;})
            << std::endl;
  //for (auto const& s: user_data.species->idx2name) {
  //  if (user_data.species->abundances[s.first] > 1.0e-40) {
  //    std::cout << s.second << " " << user_data.species->abundances[s.first] << std::endl;
  //  }
  //}

  LOGIS::loadSpeciesEnthalpies(*user_data.species, "Species_enthalpy.dat");
  std::cout << "Number of species with enthalpies: "
            << user_data.species->enthalpies.size() << std::endl;
  //for (auto const& s: user_data.species->enthalpies) {
  //  std::cout << user_data.species->idx2name[s.first] << " " << s.second << std::endl;
  //}

  CALC_RATE::assignReactionHandlers(user_data);

  TYPES::DTP_Y y = N_VClone(user_data.species->y);
  RATE_EQ::Updater_RE updater_re;
  updater_re.set_user_data(&user_data);
  updater_re.initialize_solver(0.0, user_data.species->y);
  updater_re.update(0.0, 1e7, y);

//  updater_re.
//    set_reactions().
//    set_physical_params().
//    set_chemical_params().
//    set_initial_y().
//    initialize_solver();

  return 0;
}
