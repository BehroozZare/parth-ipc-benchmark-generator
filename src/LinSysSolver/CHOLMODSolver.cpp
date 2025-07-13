//
//  CHOLMODSolver.cpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//

#ifdef IPC_WITH_CHOLMOD

#include "CHOLMODSolver.hpp"
#include "getRSS.hpp"
#include "omp.h"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace IPC {

///---------------------------------------------------------------------------------------\n
/// print_csc - save the csc sparse matrix in matrix market format (The matrix
/// is symmetric and lower triangular)
///---------------------------------------------------------------------------------------\n
void print_csc(
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
CHOLMODSolver<vectorTypeI, vectorTypeS>::CHOLMODSolver(
    std::vector<double> config_par) {
  Base::initConfig(config_par);
  cholmod_start(&cm);
  A = NULL;
  L = NULL;
  b = NULL;

  x_solve = NULL;

  x_cd = y_cd = NULL;
  Ai = Ap = Ax = NULL;
  bx = NULL;
  solutionx = x_cdx = y_cdx = NULL;
}

template <typename vectorTypeI, typename vectorTypeS>
CHOLMODSolver<vectorTypeI, vectorTypeS>::~CHOLMODSolver(void) {
  if (A) {
    A->i = Ai;
    A->p = Ap;
    A->x = Ax;
    cholmod_free_sparse(&A, &cm);
  }

  cholmod_free_factor(&L, &cm);

  if (b) {
    b->x = bx;
    cholmod_free_dense(&b, &cm);
  }

  if (x_cd) {
    x_cd->x = x_cdx;
    cholmod_free_dense(&x_cd, &cm);
  }

  if (y_cd) {
    y_cd->x = y_cdx;
    cholmod_free_dense(&y_cd, &cm);
  }

  if (x_solve) {
    cholmod_free_dense(&x_solve, &cm);
  }

  cholmod_finish(&cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::set_pattern(
    const std::vector<std::set<int>> &vNeighbor,
    const std::set<int> &fixedVert) {
  Base::set_pattern(vNeighbor, fixedVert);

  Base::load_time = omp_get_wtime();
  // TODO: directly save into A
  if (!A) {
    A = cholmod_allocate_sparse(Base::numRows, Base::numRows, Base::ja.size(),
                                true, true, -1, CHOLMOD_REAL, &cm);
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    // -1: upper right part will be ignored during computation
  }
  Base::ia.array() -= 1;
  Base::ja.array() -= 1; // CHOLMOD's index starts from 0
  A->i = Base::ja.data();
  A->p = Base::ia.data();
  A->x = Base::a.data();

  Base::parth.setMeshPointers(Base::M_n, Base::Mp.data(), Base::Mi.data());
  Base::N = Base::numRows;
  Base::NNZ = Base::a.size();
  Base::load_time = omp_get_wtime() - Base::load_time;
}
template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::set_pattern(
    const Eigen::SparseMatrix<double> &mtr) {
  Base::set_pattern(mtr);
  Base::load_time = omp_get_wtime();
  if (!A) {
    A = cholmod_allocate_sparse(Base::numRows, Base::numRows, mtr.nonZeros(),
                                true, true, -1, CHOLMOD_REAL, &cm);
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    // -1: upper right part will be ignored during computation
    A->i = Base::ja.data();
    A->p = Base::ia.data();
    A->x = Base::a.data();
  }

  Base::parth.setMeshPointers(Base::M_n, Base::Mp.data(), Base::Mi.data());
  Base::N = Base::numRows;
  Base::NNZ = Base::a.size();

  Base::load_time = omp_get_wtime() - Base::load_time;
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::load(const char *filePath,
                                                   Eigen::VectorXd &rhs) {
  Base::load(filePath, rhs);

  // TODO: directly save into A
  Base::load_time = omp_get_wtime();
  if (!A) {
    A = cholmod_allocate_sparse(Base::numRows, Base::numRows, Base::ja.size(),
                                true, true, -1, CHOLMOD_REAL, &cm);
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    // -1: upper right part will be ignored during computation
  }
  Base::ia.array() -= 1;
  Base::ja.array() -= 1; // CHOLMOD's index starts from 0
  A->i = Base::ja.data();
  A->p = Base::ia.data();
  A->x = Base::a.data();
  Base::load_time = omp_get_wtime() - Base::load_time;
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::innerAnalyze_pattern(void) {
  std::cout << "++++++++++++++++++ Cholmod: Analysing *********************"
            << std::endl;

  cholmod_free_factor(&L, &cm);

  if (Base::active_parth) {
    Base::parth_time = omp_get_wtime();
    Base::parth.computePermutation(Base::perm, Base::dim);
    Base::parth_time = omp_get_wtime() - Base::parth_time;
    assert(Base::perm.size() == Base::N);
  } else {
    Base::parth.integrator.reuse_ratio = 0;
  }

  cm.supernodal = CHOLMOD_SUPERNODAL;
  if (Base::active_parth) {
    cm.nmethods = 1;
    cm.method[0].ordering = CHOLMOD_GIVEN;
    L = cholmod_analyze_p(A, Base::perm.data(), NULL, 0, &cm);
  } else {
    cm.nmethods = 1;
    if (Base::fill_reducing_type == FillReducingType::METIS) {
      cm.method[0].ordering = CHOLMOD_METIS;
      std::cout << "*** CHOLMOD: Choosing METIS" << std::endl;
    } else if (Base::fill_reducing_type == FillReducingType::AMD) {
      cm.method[0].ordering = CHOLMOD_AMD;
      std::cout << "*** CHOLMOD: Choosing AMD" << std::endl;
    } else {
      std::cerr << "*** CHOLMOD: UNKNOWN REORDERING TYPE -> Choosing METIS"
                << std::endl;
    };
    L = cholmod_analyze(A, &cm);
    //    // Storing the permutation array
    //    int* perm = (int*)L->Perm;
    //    curr_perm.resize(Base::numRows);
    //    std::copy(perm, perm + Base::numRows, curr_perm.data());
    //    // Do some analysis with curr_perm
    //    if (prev_perm.size() != 0) {
    //        double dot_product = 0;
    //        assert(curr_perm.size() == prev_perm.size());
    //        for (int iter = 0; iter < curr_perm.size(); iter++) {
    //            dot_product += (curr_perm[iter] * prev_perm[iter]);
    //        }
    //
    //        double curr_per_norm2 = 0;
    //        for (int iter = 0; iter < curr_perm.size(); iter++) {
    //            curr_per_norm2 += (curr_perm[iter] * curr_perm[iter]);
    //        }
    //        curr_per_norm2 = std::sqrt(curr_per_norm2);
    //
    //        double prev_per_norm2 = 0;
    //        for (int iter = 0; iter < prev_perm.size(); iter++) {
    //            prev_per_norm2 += (prev_perm[iter] * prev_perm[iter]);
    //        }
    //        prev_per_norm2 = std::sqrt(prev_per_norm2);
    //
    //        double cos = dot_product / (curr_per_norm2 * prev_per_norm2);
    //        std::cout << "****** Difference " << cos << std::endl;
    //    }
    //    prev_perm = curr_perm;
  }
  Base::L_NNZ = cm.lnz * 2 - Base::numRows;
}

template <typename vectorTypeI, typename vectorTypeS>
bool CHOLMODSolver<vectorTypeI, vectorTypeS>::innerFactorize(void) {
  if (!Base::active_parth) {
    Base::saveHessian("last");
  }
  std::cout << "++++++++++++++++++ Cholmod: Factorizing *********************"
            << std::endl;
  cholmod_factorize(A, L, &cm);
  return cm.status != CHOLMOD_NOT_POSDEF;
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::innerSolve(
    Eigen::VectorXd &rhs, Eigen::VectorXd &result) {
  std::cout << "++++++++++++++++++ Cholmod: Solve *********************"
            << std::endl;
  // TODO: directly point to rhs?
  if (!b) {
    b = cholmod_allocate_dense(Base::numRows, 1, Base::numRows, CHOLMOD_REAL,
                               &cm);
    bx = b->x;
  }
  b->x = rhs.data();

  if (x_solve) {
    cholmod_free_dense(&x_solve, &cm);
  }

  x_solve = cholmod_solve(CHOLMOD_A, L, b, &cm);

  //  Base::residual_comp_time = omp_get_wtime();
  //  double one[2] = {1, 0}, m1[2] = {-1, 0};
  //  cholmod_dense *r = cholmod_copy_dense(b, &cm);
  //  cholmod_sdmult(A, 0, m1, one, x_solve, r, &cm);
  //  Base::residual = cholmod_norm_dense(r, 0, &cm);
  //  cholmod_free_dense(&r, &cm);
  //  Base::residual_comp_time = omp_get_wtime() - Base::residual_comp_time;

  result.conservativeResize(rhs.size());
  memcpy(result.data(), x_solve->x, result.size() * sizeof(result[0]));

  // Save files for testing TODO:Delete this when development is finished

  //  const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
  //                                         Eigen::DontAlignCols, ", ", "\n");

  //  std::ofstream file1(Base::outputFolderPath + "/rhs_" +
  //                      std::to_string(Base::frame) + "_" +
  //                      std::to_string(Base::iteration));
  //  if (file1.is_open()) {
  //    file1 << rhs.format(CSVFormat);
  //    file1.close();
  //  }
  //
  //  // Save files for testing TODO:Delete this when development is finished
  //  std::ofstream file2(Base::outputFolderPath + "/sol_" +
  //                      std::to_string(Base::frame) + "_" +
  //                      std::to_string(Base::iteration));
  //  if (file2.is_open()) {
  //    file2 << result.format(CSVFormat);
  //    file2.close();
  //  }
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::multiply(const Eigen::VectorXd &x,
                                                       Eigen::VectorXd &Ax) {
  assert(x.size() == Base::numRows);

  if (!x_cd) {
    x_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows, CHOLMOD_REAL,
                                  &cm);
    x_cdx = x_cd->x;
  }
  x_cd->x = (void *)x.data();

  Ax.conservativeResize(Base::numRows);
  if (!y_cd) {
    y_cd = cholmod_allocate_dense(Base::numRows, 1, Base::numRows, CHOLMOD_REAL,
                                  &cm);
    y_cdx = y_cd->x;
  }
  y_cd->x = (void *)Ax.data();

  double alpha[2] = {1.0, 1.0}, beta[2] = {0.0, 0.0};

  cholmod_sdmult(A, 0, alpha, beta, x_cd, y_cd, &cm);
}

template <typename vectorTypeI, typename vectorTypeS>
void CHOLMODSolver<vectorTypeI, vectorTypeS>::outputFactorization(
    const std::string &filePath) {
  cholmod_sparse *spm = cholmod_factor_to_sparse(L, &cm);

  FILE *out = fopen(filePath.c_str(), "w");
  assert(out);

  cholmod_write_sparse(out, spm, NULL, "", &cm);

  fclose(out);
}

template class CHOLMODSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC

#endif
