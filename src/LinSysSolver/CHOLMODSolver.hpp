//
//  CHOLMODSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//
#pragma once

#ifdef IPC_WITH_CHOLMOD
#include "LinSysSolver.hpp"
#include "Parth.h"

#include "cholmod.h"

#include <Eigen/Eigen>

#include <set>
#include <vector>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
class CHOLMODSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
  typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
  cholmod_common cm;
  cholmod_sparse *A;
  cholmod_factor *L;
  cholmod_dense *b, *solution;
  cholmod_dense *x_cd, *y_cd; // for multiply

  cholmod_dense *x_solve;

  void *Ai, *Ap, *Ax, *bx, *solutionx, *x_cdx, *y_cdx;

public:
  CHOLMODSolver(std::vector<double> config_par);
  ~CHOLMODSolver(void);

  LinSysSolverType type() const override { return LinSysSolverType::CHOLMOD; }

  void set_pattern(const std::vector<std::set<int>> &vNeighbor,
                   const std::set<int> &fixedVert) override;
  void set_pattern(
      const Eigen::SparseMatrix<double> &mtr) override; // NOTE: mtr must be SPD
  void load(const char *filePath, Eigen::VectorXd &rhs) override;

  void innerAnalyze_pattern(void) override;

  bool innerFactorize(void) override;

  void innerSolve(Eigen::VectorXd &rhs, Eigen::VectorXd &result) override;

  void multiply(const Eigen::VectorXd &x, Eigen::VectorXd &Ax) override;

  void outputFactorization(const std::string &filePath) override;
};

} // namespace IPC

#endif
