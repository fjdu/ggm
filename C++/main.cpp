#include <iostream>
#include <string>
#include <algorithm>
#include "logistics.hpp"
#include "types.hpp"
#include "rate_equation_lsode.hpp"
#include "constants.hpp"
#include "calculate_reaction_rate.hpp"

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

  LOGIS::load_reactions("reactions_test.dat", &user_data);
  std::cout << "Number of reactions: "
            << user_data.reactions->size() << std::endl;
  std::cout << "Number of species: "
            << user_data.species->name2idx.size() << std::endl;
  std::cout << "Number of reaction types: "
            << user_data.reaction_types->size() << std::endl;
  for (auto const& i: *(user_data.reaction_types)) {
    std::cout << i.first << ": " << i.second << "\n";
  }
  LOGIS::assort_reactions(*(user_data.reactions), user_data.other_data);
  std::cout << "Number of adsorption reactions: "
            << user_data.other_data->ads_reactions.size() << std::endl;
  std::cout << "Number of evaporation reactions: "
            << user_data.other_data->eva_reactions.size() << std::endl;

  LOGIS::assignElementsToSpecies(*(user_data.species), CONST::element_masses);
  LOGIS::calculateSpeciesMasses(*(user_data.species), CONST::element_masses);
  LOGIS::calculateSpeciesVibFreqs(*(user_data.species), *(user_data.reactions));
  LOGIS::calculateSpeciesDiffBarriers(*(user_data.species), *(user_data.reactions));
  LOGIS::calculateSpeciesQuantumMobilities(*(user_data.species), *(user_data.reactions));
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
  for (auto const& s: user_data.species->idx2name) {
    if (user_data.species->abundances[s.first] > 1.0e-40) {
      std::cout << s.second << " " << user_data.species->abundances[s.first] << std::endl;
    }
  }

  LOGIS::loadSpeciesEnthalpies(*user_data.species, "Species_enthalpy.dat");
  std::cout << "Number of species with enthalpies: "
            << user_data.species->enthalpies.size() << std::endl;
  //for (auto const& s: user_data.species->enthalpies) {
  //  std::cout << user_data.species->idx2name[s.first] << " " << s.second << std::endl;
  //}

  CALC_RATE::assignReactionHandlers(user_data);
  for (auto const& i: *(user_data.reaction_types)) {
    if (user_data.rate_calculators->find(i.first) ==
        user_data.rate_calculators->end()) {
      std::cout << "Reaction type " << i.first
                << " has no handler!" << std::endl;
    }
  }

  RATE_EQ::Updater_RE updater_re(user_data.species->idx2name.size());
  updater_re.set_user_data(&user_data);
  updater_re.initialize_solver(1e-6, 1e-30, 21);

  TYPES::DTP_FLOAT t = 0.0, dt=1e-10, t_ratio=1.2;
  double *y = new double[updater_re.NEQ];
  for (int i=0; i<updater_re.NEQ; ++i) {
    y[i] = user_data.species->abundances[i];
  }
  for (int i=0; i<300; ++i) {
    t = updater_re.update(t, dt, y);
    dt *= t_ratio;
    if (updater_re.ISTATE != 2) {
      std::cout << "Failed: " << updater_re.ISTATE << std::endl;
      break;
    }
  }
  std::cout << "Number of f evaluations: "
            << updater_re.IWORK[11] << std::endl;
  std::cout << "Number of jac evaluations: "
            << updater_re.IWORK[12] << std::endl;
  std::cout << "RWORK size required: "
            << updater_re.IWORK[16] << std::endl;
  std::cout << "IWORK size required: "
            << updater_re.IWORK[17] << std::endl;
  std::cout << "NNZ: "
            << updater_re.IWORK[18] << std::endl;
  //int IPIAN = updater_re.IWORK[22], IPJAN = updater_re.IWORK[23];
  //double *IWK = &(updater_re.RWORK[20]);
  //std::cout << IPIAN << " I,J " << IPJAN << std::endl;
  //std::cout << IWK[IPIAN-1] << std::endl;
  //for (int i=0; i<=updater_re.NEQ; ++i) {
  //  std::cout << i+1 << " IWK_I " << IWK[IPIAN + i] << std::endl;
  //}
  //for (int i=0; i<updater_re.NNZ; ++i) {
  //  std::cout << i+1 << " IWK_J " << IWK[IPJAN + i] << std::endl;
  //}
  //for (int i=30; i<updater_re.LIW; ++i) {
  //  std::cout << i+1 << " " << updater_re.IWORK[i] << std::endl;
  //}
  //for (int i=0; i<user_data.species->name2idx.size(); ++i) {
  //  std::cout << user_data.species->idx2name[i] << " "
  //            << user_data.species->abundances[i] << " "
  //            << user_data.species->abundances[i] - user_data.species->abundances[i]
  //            << std::endl;
  //}
  //double *ydot = new double[updater_re.NEQ];
  //double *ydot1 = new double[updater_re.NEQ];
  //for (int i=0; i<updater_re.NEQ; ++i) {
  //  double dy = y[i] * 1e-4 + 1e-30;
  //  for (int j=0; j<updater_re.NEQ; ++j) {
  //    ydot[j] = 0.0; ydot1[j] = 0.0;
  //  }
  //  updater_re.f(&updater_re.NEQ, &t, y, ydot);
  //  y[i] += dy;
  //  updater_re.f(&updater_re.NEQ, &t, y, ydot1);
  //  y[i] -= dy;
  //  for (int j=0; j<updater_re.NEQ; ++j) {
  //    double jc = (ydot1[j] - ydot[j]) / dy;
  //    std::cout << "r" << j << " c" << i << " " << jc << std::endl;
  //  }
  //}

//  updater_re.
//    set_reactions().
//    set_physical_params().
//    set_chemical_params().
//    set_initial_y().
//    initialize_solver();

  return 0;
}
