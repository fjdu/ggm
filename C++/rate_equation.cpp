#include <iostream>
#include "types.hpp"
#include "update.hpp"
#include "cvodes/cvodes.h"
#include <nvector/nvector_serial.h>         /* access to serial N_Vector            */
#include <sundials/sundials_types.h>        /* defs. of realtype, sunindextype      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvodes/cvodes_direct.h>      /* access to CVDls interface            */

#ifndef RATE_EQ_H
#define RATE_EQ_H

#define Ith(v,i)    NV_Ith_S(v,i)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i,j) /* IJth numbers rows,cols 1..NEQ */


namespace RATE_EQ {

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  TYPES::User_data *data;
  data = (TYPES::User_data*) user_data;

  for (auto const& reaction: *data->reactions) {
    if ((reaction.itype == 67)) {
      continue;
    }
    //std::cout << "reaction.itype: " << reaction.itype << std::endl;
    //std::cout << "t = " << t << std::endl;
    //std::cout << "y[0] = " << NV_Ith_S(y, 0) << std::endl;
    TYPES::DTP_FLOAT r =
      ((*(data->rate_calculators))[reaction.itype])(t, y,
        reaction, *(data->physical_params),
        *(data->species), *(data->other_data));
    //for (auto const& i: reaction.idxReactants) {
    //  std::cout << data->species->idx2name[i] << " ";
    //}
    //std::cout <<  " -> ";
    //for (auto const& i: reaction.idxProducts) {
    //  std::cout << data->species->idx2name[i] << " ";
    //}
    //std::cout <<  " " << r << std::endl;
    for (auto const& i: reaction.idxReactants) {
      NV_Ith_S(ydot, i) -= r;
    }
    for (auto const& i: reaction.idxProducts) {
      NV_Ith_S(ydot, i) += r;
    }
  }

  return 0;
}


//static int jtv(N_Vector v, N_Vector Jv, realtype t,
//               N_Vector y, N_Vector fy,
//               void *user_data, N_Vector tmp)
//{
//  TYPES::User_data *data;
//  data = (TYPES::User_data*) user_data;
//
//  for (int i=0; i<NV_LENGTH_S(Jv); ++i) {
//    NV_Ith_S(Jv, i) = 0.0;
//  }
//
//  for (auto const& reaction: *data->reactions) {
//    if ((reaction.itype == 67)) {
//      continue;
//    }
//    std::vector<TYPES::DTP_FLOAT> drdy = ((*(data->drdy_calculators))[reaction.itype])(
//      t, y, reaction, *(data->physical_params), *(data->species), *(data->other_data));
//
//    TYPES::DTP_FLOAT jdotv = 0.0;
//    for (int i=0; i<reaction.idxReactants.size(); ++i) {
//      jdotv += drdy[i] * NV_Ith_S(v, reaction.idxReactants[i]);
//    }
//
//    for (auto const& i: reaction.idxReactants) {
//      NV_Ith_S(Jv, i) -= jdotv;
//    }
//    for (auto const& i: reaction.idxProducts) {
//      NV_Ith_S(Jv, i) += jdotv;
//    }
//  }
//
//  return 0;
//}


static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  TYPES::User_data *data;
  data = (TYPES::User_data*) user_data;

  SUNMatZero(J);

  for (auto const& r: *(data->reactions)) {
    if (r.itype == 67) {
      continue;
    }
    std::vector<TYPES::DTP_FLOAT> drdy =
      ((*(data->drdy_calculators))[r.itype])(t, y, r,
        *(data->physical_params), *(data->species), *(data->other_data));
    for (int i=0; i < r.idxReactants.size(); ++i) {
      for (auto const& j: r.idxReactants) {
        IJth(J, j, r.idxReactants[i]) -= drdy[i];
      }
      for (auto const& j: r.idxProducts) {
        IJth(J, j, r.idxReactants[i]) += drdy[i];
      }
    }
  }
  //for (int i=0; i<NV_LENGTH_S(y); ++i) {
  //  for (int j=0; j<NV_LENGTH_S(y); ++j) {
  //    std::cout << IJth(J, i, j) << " ";
  //  }
  //  std::cout << std::endl;
  //}

  return(0);
}


class Updater_RE: Update::Updater {
  public:
    Updater_RE() {
      solver_initialized = false;
    }

    TYPES::DTP_FLOAT update(const TYPES::DTP_FLOAT t, const TYPES::DTP_FLOAT dt, TYPES::DTP_Y y);

    int set_user_data(TYPES::User_data *_u);

    Updater_RE& set_recorder();
    Updater_RE& initialize_solver(TYPES::DTP_FLOAT t0, TYPES::DTP_Y y0,
        TYPES::DTP_FLOAT reltol, TYPES::DTP_FLOAT absto, int lmm);

    ~Updater_RE() {
      if (solver_initialized) {
        SUNLinSolFree(LS);
        CVodeFree(&cvode_mem);
      }
      if (user_data != NULL) {
        //std::cout << user_data->species->massSpecies.size() << std::endl;
      }
    }
  private:
    bool solver_initialized;
    void *cvode_mem;
    TYPES::User_data *user_data;
    SUNLinearSolver LS;
    TYPES::DTP_Y y0;
};


TYPES::DTP_FLOAT Updater_RE::update(const TYPES::DTP_FLOAT t,
    const TYPES::DTP_FLOAT dt,
    TYPES::DTP_Y y)
{
  int itask = CV_NORMAL;
  TYPES::DTP_FLOAT tret;
  int retval = CVode(cvode_mem, t+dt, y, &tret, itask);

  if (retval == CV_SUCCESS) {
  }
  std::cout << "t reaches: " << tret << std::endl;

  return tret;
}


int Updater_RE::set_user_data(TYPES::User_data *_u) {
  user_data = _u;
  return 0;
}


Updater_RE& Updater_RE::initialize_solver(
    TYPES::DTP_FLOAT t0,
    TYPES::DTP_Y y0,
    TYPES::DTP_FLOAT reltol=1e-6,
    TYPES::DTP_FLOAT abstol=1e-30,
    int lmm=CV_BDF)
{
  sunindextype NEQ = NV_LENGTH_S(y0);
  SUNMatrix A;
  N_Vector constraints = N_VNew_Serial(NEQ);
  N_VConst(1.0, constraints);

  cvode_mem = CVodeCreate(lmm);
  if (cvode_mem == NULL) {
    std::cout << "Error calling CVodeCreate." << std::endl;
    return *this;
  }

  int flag = CVodeInit(cvode_mem, f, t0, y0);
  if (flag != CV_SUCCESS) {
    std::cout << "Error calling CVodeInit: " << flag << std::endl;
    return *this;
  }

  flag = CVodeSetUserData(cvode_mem, user_data);
  if (flag != CV_SUCCESS) {
    std::cout << "Error calling CVodeSetUserData: " << flag << std::endl;
    return *this;
  }

  flag = CVodeSStolerances(cvode_mem, reltol, abstol);

  A = SUNDenseMatrix(NEQ, NEQ);
  LS = SUNLinSol_Dense(y0, A);

  flag = CVodeSetLinearSolver(cvode_mem, LS, A);
  flag = CVodeSetJacFn(cvode_mem, Jac);
  CVodeSetConstraints(cvode_mem, constraints);
  N_VDestroy(constraints);

  solver_initialized = true;
  return *this;
}


}

#endif //RATE_EQ_H
