#include <iostream>
#include <string>
#include <algorithm>
#include "logistics.hpp"
#include "types.hpp"
#include "rate_equation_lsode.hpp"
#include "constants.hpp"
#include "calculate_reaction_rate.hpp"
#include <ctime>

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

  LOGIS::load_reactions("rate06_withgrain_lowH2Bind_hiOBind_lowCObind.dat", &user_data);
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
  for (int i=0; i<user_data.species->idx2name.size(); ++i) {
    if (user_data.species->abundances[i] > 1.0e-40) {
      std::cout << user_data.species->idx2name[i] << " "
                << user_data.species->abundances[i] << std::endl;
    }
  }

  LOGIS::loadSpeciesEnthalpies(*user_data.species, "Species_enthalpy.dat");
  std::cout << "Number of species with enthalpies: "
            << user_data.species->enthalpies.size() << std::endl;
  //for (auto const& s: user_data.species->enthalpies) {
  //  std::cout << user_data.species->idx2name[s.first]
  //            << " " << s.second << std::endl;
  //}

  CALC_RATE::assignReactionHandlers(user_data);
  for (auto const& i: *(user_data.reaction_types)) {
    if (user_data.rate_calculators->find(i.first) ==
        user_data.rate_calculators->end()) {
      std::cout << "Reaction type " << i.first
                << " has no handler!" << std::endl;
    }
  }

  TYPES::Recorder recorder("evol_001.dat");
  recorder.write_header(user_data.species->idx2name);

  RATE_EQ::Updater_RE updater_re(user_data.species->idx2name.size());
  updater_re.set_user_data(&user_data);
  updater_re.initialize_solver(1e-6, 1e-30);

  TYPES::DTP_FLOAT t = 0.0, dt=1e-2, t_ratio=1.1;
  double *y = new double[updater_re.NEQ];
  for (int i=0; i<updater_re.NEQ; ++i) {
    y[i] = user_data.species->abundances[i];
  }

  recorder.write_row(t, updater_re.NEQ, y);
  std::cout << std::endl;

  clock_t rt_begin = std::clock();
  for (int i=0; i<400; ++i) {
    t = updater_re.update(t, dt, y);
    recorder.write_row(t/CONST::phy_SecondsPerYear, updater_re.NEQ, y);
    dt *= t_ratio;
    if (updater_re.ISTATE != 2) {
      std::cout << "Failed: " << updater_re.ISTATE << std::endl;
      if ((updater_re.ISTATE == -1) ||
          (updater_re.ISTATE == -4) ||
          (updater_re.ISTATE == -5)) {
        updater_re.ISTATE = 3;
      } else {
        break;
      }
    }
  }
  clock_t rt_end = std::clock();
  double elapsed_secs = double(rt_end - rt_begin) / CLOCKS_PER_SEC;
  std::cout << elapsed_secs << " seconds elapsed." << std::endl;

  delete [] y;

  return 0;
}
