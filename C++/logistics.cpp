#ifndef LOGIS_H
#define LOGIS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <regex>
#include <set>
#include "types.hpp"
#include "utils.cpp"
#include "math.h"

void load_enthalpies();
void load_initial_condition();
void load_config();
void record();
void solve();

namespace LOGIS {


int load_reactions(const std::string& fname,
    TYPES::User_data *user_data,
    int nReactants=3, int nProducts=4, int nABC=3,
    int lenSpeciesName=12, int lenABC=9, int nT=2, int lenT=6, int lenType=3,
    int rowlen_min=126)
{
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");

  TYPES::Species& species = *user_data->species;
  TYPES::Reactions& reactions = *user_data->reactions;
  std::set<int>& reaction_types = *user_data->reaction_types;

  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {
      if (std::regex_match(line, comment)) {
        // std::cout << "Comments!" << std::endl;
        continue;
      }
      if (std::regex_match(line, emptyline)) {
        // std::cout << "Empty line!" << std::endl;
        continue;
      }
      if (line.size() < rowlen_min) {
        // std::cout << "Too short line!" << std::endl;
        continue;
      }

      TYPES::Reaction reaction;

      for (int i=0; i<nReactants+nProducts; ++i) {
        std::string name = line.substr(i*lenSpeciesName, lenSpeciesName);
        name = UTILS::trim(name);

        if ((name.size() > 0) &&
            (species.name2idx.find(name) == species.name2idx.end())) {
          int nSpecies = species.name2idx.size();
          species.name2idx[name] = nSpecies;
          species.idx2name[nSpecies] = name;
        }

        if (name.size() > 0) {
          int iSpecies = species.name2idx[name];
          if (i < nReactants) {
            reaction.idxReactants.push_back(iSpecies);
          } else {
            reaction.idxProducts.push_back(iSpecies);
          }
        }
      }

      int iABC = (nReactants + nProducts) * lenSpeciesName;
      for (int i=0; i<nABC; ++i) {
        std::string tmp = line.substr(iABC+i*lenABC, lenABC);
        tmp = UTILS::trim(tmp);
        if (tmp.size() > 0) {
          reaction.abc.push_back(std::stod(tmp));
        } else {
          reaction.abc.push_back(0.0);
        }
      }

      int iT = iABC + nABC * lenABC;
      for (int i=0; i<nT; ++i) {
        std::string tmp = line.substr(iT+i*lenT, lenT);
        tmp = UTILS::trim(tmp);
        if (tmp.size() > 0) {
          reaction.Trange.push_back(std::stod(tmp));
        } else {
          reaction.Trange.push_back(0.0);
        }
      }

      int iType = iT + nT * lenT;
      std::string tmp = line.substr(iType, lenType);
      reaction.itype = std::stoi(tmp);
      reaction_types.insert(reaction.itype);

      reactions.push_back(reaction);

      //for (auto const& i: reaction.idxReactants) {
      //  std::cout << species.idx2name[i] << " ("<< i << ") ";
      //}
      //std::cout << " -> ";
      //for (auto const& i: reaction.idxProducts) {
      //  std::cout << species.idx2name[i] << " (" << i << ") ";
      //}
      //std::cout << " " << reaction.abc[0] << " " << reaction.abc[1] << " " << reaction.abc[2]
      //          << " " << reaction.itype
      //          << " " << reaction.Trange[0] << " " << reaction.Trange[1]
      //          << std::endl;
    }
  } else {
    return -1;
  }

  inputFile.close();

  return 0;
}


std::map<std::string, int> assignElementsToOneSpecies(
    const std::string& name, const TYPES::Elements& elements)
{
  std::regex npat(R"(^(\d+))");

  std::map<std::string, int> eleDict;
  for (auto const& e: elements) {
    eleDict[e.first] = 0;
  }

  for (int iBg=0; iBg<name.size();) {
    bool found = false;
    for (int nlen=name.size()-iBg; nlen>0; --nlen) {
      auto q = elements.find(name.substr(iBg, nlen));
      if (q == elements.end()) {
        continue;
      }

      iBg += q->first.size();

      std::smatch matches;
      if (std::regex_search(name.substr(iBg), matches, npat)) {
        eleDict[q->first] += std::stoi(matches.str(1));
        iBg += matches.str(1).size();
      } else {
        eleDict[q->first] += 1;
      }
      found = true;
      break;
    }
    if (!found) {
      ++iBg;
    }
  }
  return eleDict;
}


int assignElementsToSpecies(TYPES::Species& species, TYPES::Elements& elements) {
  for (auto const& s: species.name2idx) {
    species.elementsSpecies[s.second] = assignElementsToOneSpecies(s.first, elements);
  }
  return 0;
}


int classifySpeciesByPhase(TYPES::Species& species) {
  for (auto const& s: species.name2idx) {
    if (s.first[0] == 'm') {
      species.mantleSpecies.insert(s.first);
    } else if (s.first[0] == 'g') {
      species.surfaceSpecies.insert(s.first);
    } else {
      species.gasSpecies.insert(s.first);
    }
  }
  return 0;
}


int calculateSpeciesMasses(TYPES::Species& species, TYPES::Elements& elements) {
  for (auto const& s: species.name2idx) {
    species.massSpecies[s.second] = 0.0;
    auto const& eleDict = species.elementsSpecies[s.second];
    for (auto const& e: eleDict) {
      species.massSpecies[s.second] += elements[e.first] * ((double)e.second);
    }
  }
  return 0;
}


int calculateSpeciesVibFreqs(TYPES::Species& species, TYPES::Reactions& reactions) {
  for (auto const& r: reactions) {
    if (r.itype == 62) {
      species.vibFreqs[r.idxReactants[0]] =
        sqrt(2.0 * CONST::phy_SitesDensity_CGS
           * CONST::phy_kBoltzmann_CGS * r.abc[2]
           / (CONST::PI * CONST::PI)
           / (CONST::phy_mProton_CGS * species.massSpecies[r.idxReactants[0]]));
    }
  }
  return 0;
}


int loadInitialAbundances(TYPES::Species& species, std::string fname) {
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s+(\S+))");

  if (species.abundances.size() != species.name2idx.size()) {
    species.abundances = std::vector<TYPES::DTP_FLOAT>(species.name2idx.size());
  }
  for (int i=0; i<species.abundances.size(); ++i) {
    species.abundances[i] = 0.0;
  }

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
        std::string name = matches.str(1);
        if (species.name2idx.find(name) != species.name2idx.end()) {
          TYPES::DTP_FLOAT val = std::stod(matches.str(2));
          species.abundances[species.name2idx[name]] = val;
        } else {
          std::cout << "Unrecognized line: " << line << std::endl;
        }
      }
    }
  } else {
    std::cout << "Error reading file: " << fname << std::endl;
  }
  species.y = N_VNew_Serial(species.abundances.size());
  for (int i=0; i<species.abundances.size(); ++i) {
    NV_Ith_S(species.y, i) = species.abundances[i];
  }
  return 0;
}


int loadSpeciesEnthalpies(TYPES::Species& species, std::string fname) {
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s+(\S+))");

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
        std::string name = matches.str(1);
        if (species.name2idx.find(name) != species.name2idx.end()) {
          TYPES::DTP_FLOAT val = std::stod(matches.str(2));
          species.enthalpies[species.name2idx[name]] = val;
        } else {
          //std::cout << "Invalid line: " << line << std::endl;
        }
      }
    }
  } else {
    std::cout << "Error reading file: " << fname << std::endl;
  }
  return 0;
}


}

#endif // LOGIS_H
