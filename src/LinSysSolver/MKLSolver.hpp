//
//  CHOLMODSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//
#pragma once

#ifndef IPC_MKL_SOLVER
#define IPC_MKL_SOLVER

#include "Parth.h"

#ifdef IPC_WITH_MKL

#include "LinSysSolver.hpp"
#include "mkl.h"

#include <Eigen/Eigen>

#include <set>
#include <vector>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
class MKLSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
  typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
  MKL_INT mtype = 2; /* Real unsymmetric matrix */
  MKL_INT nrhs = 1;  /* Number of right hand sides. */
  long int pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  MKL_INT maxfct, mnum, phase, error, msglvl;
  /* Auxiliary variables. */
  double ddum;  /* Double dummy */
  MKL_INT idum; /* Integer dummy. */
  Eigen::SparseMatrix<double, Eigen::ColMajor> A_CSC;
  MKL_INT *Ap;
  MKL_INT *Ai;
  double *Ax;
  MKL_INT N_MKL;

public:
  MKLSolver(std::vector<double> config_par);
  ~MKLSolver(void);

  LinSysSolverType type() const override { return LinSysSolverType::MKL; }

  void set_pattern(const std::vector<std::set<int>> &vNeighbor,
                   const std::set<int> &fixedVert) override;
  void set_pattern(
      const Eigen::SparseMatrix<double> &mtr) override; // NOTE: mtr must be SPD

  void innerAnalyze_pattern(void) override;

  bool innerFactorize(void) override;

  void innerSolve(Eigen::VectorXd &rhs, Eigen::VectorXd &result) override;
};

} // namespace IPC

#endif
#endif