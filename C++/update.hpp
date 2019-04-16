// Given (t, y, physical_params, chem_params, reactions, updater), retrieve y(t+dt).

#ifndef UPDATE_H
#define UPDATE_H

#include <string>
#include "types.hpp"

namespace Update {

class Updater {
  public:
    virtual TYPES::DTP_FLOAT update(
        const TYPES::DTP_FLOAT t,
        const TYPES::DTP_FLOAT dt,
        TYPES::DTP_Y y) = 0;
};


}
#endif //UPDATE_H
