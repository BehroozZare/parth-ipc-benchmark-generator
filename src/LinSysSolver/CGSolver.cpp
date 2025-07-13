//
//  EigenLibSolver.cpp
//  IPC
//
//  Created by Minchen Li on 6/30/18.
//

#include "CGSolver.hpp"
#include "omp.h"
#include <fstream>
#include <iomanip> // std::setprecision

namespace IPC {

void print_cscCG(
    size_t n,                                ///<[in] Number of nodes
    int *Ap,                                 ///<[in] pointer array
    int *Ai,                                 ///<[in] index array
    double *Ax,                              ///<[in] value array
    std::streambuf *out = std::cout.rdbuf(), ///<[in] output buffer
    const std::string indent = "  ",         ///<[in] indent definition
    const std::string &beg =
        "%%MatrixMarket matrix coordinate real symmetric" ///<[in] Header of the
) {
  if (Ap == NULL || Ai == NULL) {
    std::cout << "The matrix is not defined" << std::endl;
    return;
  }
  std::streambuf *sb_cout_backup = std::cout.rdbuf();
  std::cout.rdbuf(out);
  std::cout << indent << beg << "\n";
  size_t nnz = n > 0 ? Ap[n] : 0;
  std::cout << indent << n << " " << n << " " << nnz << "\n";
  for (auto i = 0; i < n; ++i) {
    for (auto j = Ap[i]; j < Ap[i + 1]; ++j) {
      assert(j < nnz);
      if (Ax == NULL) {
        std::cout << indent << Ai[j] + 1 << " " << i + 1 << " "
                  << std::setprecision(20) << 1;
      } else {
        std::cout << indent << Ai[j] + 1 << " " << i + 1 << " "
                  << std::setprecision(20) << Ax[j];
      }

      std::cout << "\n";
    }
  }
  std::cout.rdbuf(sb_cout_backup);
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::set_pattern(
    const std::vector<std::set<int>> &vNeighbor,
    const std::set<int> &fixedVert) {
  Base::set_pattern(vNeighbor, fixedVert);

  // TODO: directly save into mtr
  Base::ia.array() -= 1.0;
  Base::ja.array() -= 1.0;
  coefMtr.resize(Base::numRows, Base::numRows);
  coefMtr.reserve(Base::ja.size());

  memcpy(coefMtr.outerIndexPtr(), Base::ia.data(),
         Base::ia.size() * sizeof(Base::ia[0]));
  memcpy(coefMtr.innerIndexPtr(), Base::ja.data(),
         Base::ja.size() * sizeof(Base::ja[0]));
  memcpy(coefMtr.valuePtr(), Base::a.data(),
         Base::a.size() * sizeof(Base::a[0]));
}
template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::set_pattern(
    const Eigen::SparseMatrix<double> &mtr) // NOTE: mtr must be SPD
{
  Base::numRows = static_cast<int>(mtr.rows());
  coefMtr = mtr;
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::innerAnalyze_pattern(void) {
  Base::analyze_time = omp_get_wtime();
  CG.analyzePattern(coefMtr);
  Base::analyze_time = omp_get_wtime() - Base::analyze_time;
  Base::total_analyze_time += Base::analyze_time;
}

template <typename vectorTypeI, typename vectorTypeS>
bool CGSolver<vectorTypeI, vectorTypeS>::innerFactorize(void) {
  std::cout << "======================================= CG is being used "
               "=================================="
            << std::endl;
  Base::saveHessian("last");
  std::string uniq_f_name =
      std::to_string(Base::frame) + "_" + std::to_string(Base::iteration);
  std::string f_name =
      Base::outputFolderPath + "/mesh_" + uniq_f_name + "_IPC.mtx";
  std::ofstream f2_out(f_name);
  print_cscCG(Base::get_numCol_1D(), Base::get_p_1D(), Base::get_i_1D(),
              nullptr, f2_out.rdbuf(), "");
  f2_out.close();

  CG.setMaxIterations(1024);
  CG.setTolerance(1e-16f);
  bool succeeded = false;
  Base::factor_time = omp_get_wtime();
  CG.factorize(coefMtr);
  Base::factor_time = omp_get_wtime() - Base::factor_time;
  Base::total_factor_time += Base::factor_time;
  succeeded = (CG.info() == Eigen::Success);
  assert(succeeded);
  return succeeded;
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::innerSolve(Eigen::VectorXd &rhs,
                                                    Eigen::VectorXd &result) {
  Base::solve_time = omp_get_wtime();
  result = CG.solve(rhs);
  Base::solve_time = omp_get_wtime() - Base::solve_time;
  Base::total_solve_time += Base::solve_time;
}

template <typename vectorTypeI, typename vectorTypeS>
double CGSolver<vectorTypeI, vectorTypeS>::coeffMtr(int rowI, int colI) const {
  return Base::coeffMtr(rowI, colI);
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::setZero(void) {
  // TODO: directly manipulate valuePtr without a
  Base::setZero();
  memcpy(coefMtr.valuePtr(), Base::a.data(),
         Base::a.size() * sizeof(Base::a[0]));
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::setCoeff(int rowI, int colI,
                                                  double val) {
  // TODO: directly manipulate valuePtr without a

  if (rowI <= colI) {
    assert(rowI < Base::IJ2aI.size());
    const auto finder = Base::IJ2aI[rowI].find(colI);
    assert(finder != Base::IJ2aI[rowI].end());
    Base::a[finder->second] = val;
    coefMtr.valuePtr()[finder->second] = val;
  }
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::addCoeff(int rowI, int colI,
                                                  double val) {
  // TODO: directly manipulate valuePtr without a

  if (rowI <= colI) {
    assert(rowI < Base::IJ2aI.size());
    const auto finder = Base::IJ2aI[rowI].find(colI);
    assert(finder != Base::IJ2aI[rowI].end());
    Base::a[finder->second] += val;
    coefMtr.valuePtr()[finder->second] += val;
  }
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::setUnit_row(int rowI) {
  for (const auto &colIter : Base::IJ2aI[rowI]) {
    coefMtr.valuePtr()[colIter.second] = (colIter.first == rowI);
  }
}

template <typename vectorTypeI, typename vectorTypeS>
void CGSolver<vectorTypeI, vectorTypeS>::setUnit_col(
    int colI, const std::set<int> &rowVIs) {
  for (const auto &rowVI : rowVIs) {
    for (int dimI = 0; dimI < DIM; ++dimI) {
      int rowI = rowVI * DIM + dimI;
      if (rowI <= colI) {
        const auto finder = Base::IJ2aI[rowI].find(colI);
        if (finder != Base::IJ2aI[rowI].end()) {
          coefMtr.valuePtr()[finder->second] = (rowI == colI);
        }
      }
    }
  }
}

template class CGSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC
