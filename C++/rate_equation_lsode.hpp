#include <iostream>
#include <vector>
#include <algorithm>
#include "types.hpp"

#ifndef RATE_EQ_H
#define RATE_EQ_H

extern "C" {

void dlsodes_w(
  void (*)(int *, double *, double *, double *), //f
  int *NEQ, double *y, double *t, double *tout, int *ITOL,
  double *RTOL, double *ATOL, int *ITASK, int *ISTATE, int *IOPT,
  double *RWORK, int *LRW, int *IWORK, int *LIW,
  void (*)(int *, double *, double *, int *, double *, double *, double *),
  int *MF);

}

namespace RATE_EQ {

class Updater_RE {
  public:
    static void f(int *neq, double *t, double *y, double *ydot);
    static void jac(int *neq, double *t, double *y, int *j, double *ian, double *jan, double *pdj);

    void set_user_data(TYPES::User_data *data_) {
      data = data_;
      std::cout << "User_data set." << std::endl;
    }

    Updater_RE(int neq_) {
      NEQ = neq_;
      solver_initialized = false;
    
      for (int i=0; i<NEQ; ++i) {
        sparseMaskJac.push_back(std::vector<bool>(NEQ, false));
      }

      std::cout << "Updater_RE constructed." << std::endl;
    }

    ~Updater_RE() {
      if (solver_initialized) {
      }
      if (data != nullptr) {
      }
      if (IWORK != nullptr) {
        delete [] IWORK;
      }
      if (RWORK != nullptr) {
        delete [] RWORK;
      }
    }

    TYPES::DTP_FLOAT update(double t, double dt, double *y);

    Updater_RE& set_recorder();

    int initialize_solver(double reltol=1e-6, double abstol=1e-30, int mf=21, int LRW_fact=4);
    int makeSparse(const TYPES::Reactions& reactions,
                   std::vector<std::vector<bool> >& sps);

    static TYPES::User_data *data;
    static int *IWORK;
    static double *RWORK;
    int NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF, NNZ;
    double RTOL, ATOL;
    std::vector<std::vector<bool> > sparseMaskJac;
    bool solver_initialized;
};


}

#endif //RATE_EQ_H
