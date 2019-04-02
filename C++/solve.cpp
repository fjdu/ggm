#ifndef SOLVE_H
#define SOLVE_H

#include <map>
#include <vector>
#include <string>
#include "types.hpp"
#include "sundials/sundials_types.h"
#include "sundials/sundials_config.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_sparse.h"

namespace SOLVE {


//int solve(vector<DTP_FLOAT> y0, Recorder recorder, Analyzer analyzer, Reactions reactions, ODE_Solver ode_solver)
//{
//  DTP_FLOAT t, tout;
//  vector<DTP_FLOAT> y=y0;
//  while (recorder.get_time(&t, &tout)) {
//    y = ode_solver.solve(t, tout, y0, f, jac);
//    recorder.record(t, y);
//    analyzer.analyse(t, y, reactions);
//  }
//}




}
#endif //SOLVE_H
