//
//  CHOLMODSolver.cpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//

#ifdef IPC_WITH_MKL

#include "MKLSolver.hpp"
#include "getRSS.hpp"
#include "omp.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
MKLSolver<vectorTypeI, vectorTypeS>::MKLSolver(std::vector<double> config_par) {

  Base::initConfig(config_par);

  for (int i = 0; i < 64; i++) {
    pt[i] = 0;
    iparm[i] = 0;
  }
  //    iparm[0] = 1;
  //    iparm[1] = 2; // METIS REORDERING
  //    iparm[10] = 0; // Scaling
  //    iparm[12] = 0; // Matching
  //    iparm[17] = -1;
  //    iparm[18] = -1;
  //    iparm[20] = 0; // diagonal 1*1 pivoting
  //    iparm[26] = 1;
  //    iparm[34] = 1; // 0-based indexing
  //    iparm[36] = 0; // Blocking

  iparm[0] = 1; /* No solver default */
  if (Base::fill_reducing_type == FillReducingType::METIS) {
    iparm[1] = 2; /* Fill-in reordering from METIS */
  } else if (Base::fill_reducing_type == FillReducingType::AMD) {
    iparm[1] = 0; /* Fill-in reordering from AMD */
  }
  iparm[2] = 0;
  iparm[3] = 0; /* No iterative-direct algorithm */
  if (Base::active_parth) {
    iparm[4] = 1; /* No iterative-direct algorithm */
  } else {
    iparm[4] = 0; /* No iterative-direct algorithm */
  }
  iparm[5] = 0;  /* Write solution into x */
  iparm[6] = 0;  /* Not in use */
  iparm[7] = 1;  /* Max numbers of iterative refinement steps */
  iparm[8] = 0;  /* Not in use */
  iparm[9] = 0;  /* Perturb the pivot elements with 1E-8 */
  iparm[10] = 0; /* Use nonsymmetric permutation and scaling MPS */
  iparm[11] = 0; /* A^TX=B */
  iparm[12] = 0; /* Maximum weighted matching algorithm is switched-off (default
                    for symmetric). Try iparm[12] = 1 in case of inappropriate
                    accuracy */
  iparm[13] = 0; /* Output: Number of perturbed pivots */
  iparm[14] = 0; /* Not in use */
  iparm[15] = 0; /* Not in use */
  iparm[16] = 0; /* Not in use */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[18] = -1; /* Output: Mflops for LU factorization */
  iparm[19] = 0;  /* Output: Numbers of CG Iterations */
  iparm[20] = 1;  /*using bunch kaufman pivoting*/
  iparm[55] = 0;  /*Diagonal and pivoting control., default is zero*/
  iparm[59] =
      1; /* Use in-core intel MKL pardiso if the sze of the problem do not
            exceed MKL_PARDISO_OOC_MAX_CORE_SIZE. Otherwise, use out-of-core*/
  //
  iparm[26] = 1;
  // iparm[23] = 1; //TODO: Later enable to se if the parallelism is better
  iparm[34] = 1;
  // Because iparm[4]==0 so:
  iparm[30] = 0;
  iparm[35] = 0;

  maxfct = 1; /* Maximum number of numerical factorizations. */
  mnum = 1;   /* Which factorization to use. */
  msglvl = 0; /* Print statistical information in file */
  error = 0;  /* Initialize error flag */
  nrhs = 1;   /* Number of right hand sides. */
  mtype = 2;  // Real and SPD matrices
}

template <typename vectorTypeI, typename vectorTypeS>
MKLSolver<vectorTypeI, vectorTypeS>::~MKLSolver() {
  phase = -1; /* Release internal memory. */
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &N_MKL, &ddum, Ap, Ai, NULL,
          &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
}

template <typename vectorTypeI, typename vectorTypeS>
void MKLSolver<vectorTypeI, vectorTypeS>::set_pattern(
    const std::vector<std::set<int>> &vNeighbor,
    const std::set<int> &fixedVert) {
  Base::set_pattern(vNeighbor, fixedVert);
  Base::ia.array() -= 1;
  Base::ja.array() -= 1; // CHOLMOD's index starts from 0

  Ap = Base::ia.data();
  Ai = Base::ja.data();
  Ax = Base::a.data();
  N_MKL = Base::numRows;
  Base::N = Base::numRows;
  Base::NNZ = Base::a.size();
  Base::parth.setMeshPointers(Base::M_n, Base::Mp.data(), Base::Mi.data());
}
template <typename vectorTypeI, typename vectorTypeS>
void MKLSolver<vectorTypeI, vectorTypeS>::set_pattern(
    const Eigen::SparseMatrix<double> &mtr) {
  Base::numRows = static_cast<int>(mtr.rows());
  Base::set_pattern(mtr);
  N_MKL = Base::numRows;
  Base::N = Base::numRows;
  Base::NNZ = Base::a.size();
  Base::parth.setMeshPointers(Base::M_n, Base::Mp.data(), Base::Mi.data());
}

template <typename vectorTypeI, typename vectorTypeS>
void MKLSolver<vectorTypeI, vectorTypeS>::innerAnalyze_pattern(void) {
  std::cout << "++++++++++++++++++ MKL: Analysing *********************"
            << std::endl;
  assert(Base::N == N_MKL);
  phase = -1;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &N_MKL, Ax, Ap, Ai, NULL, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error);

  if (Base::active_parth) {
    Base::parth_time = omp_get_wtime();
    Base::parth.computePermutation(Base::perm, Base::dim);
    Base::parth_time = omp_get_wtime() - Base::parth_time;
    assert(Base::perm.size() == Base::N);
  } else {
    Base::parth.integrator.reuse_ratio = 0;
  }

  phase = 11;
  if (Base::active_parth) {
    assert(Base::perm.size() == Base::N);
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &N_MKL, Ax, Ap, Ai,
            Base::perm.data(), &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
  } else {
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &N_MKL, Ax, Ap, Ai, NULL, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
  }
  if (error != 0) {
    std::cerr << "ERROR during symbolic factorization - code: " << error
              << std::endl;
  }
}

template <typename vectorTypeI, typename vectorTypeS>
bool MKLSolver<vectorTypeI, vectorTypeS>::innerFactorize(void) {
  std::cout << "++++++++++++++++++ MKL: Factorizing *********************"
            << std::endl;
  assert(Base::N == N_MKL);
  phase = 0;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &N_MKL, Ax, Ap, Ai, NULL, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error);

  phase = 22;
  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &N_MKL, Ax, Ap, Ai, NULL, &nrhs,
          iparm, &msglvl, &ddum, &ddum, &error);

  Base::L_NNZ = iparm[17];
  if (error != 0) {
    std::cerr << "ERROR during numerical factorization - code: " << error
              << std::endl;
    return false;
  }
  return true; // TODO:CHECK FOR SPD FLAGS LATER
}

template <typename vectorTypeI, typename vectorTypeS>
void MKLSolver<vectorTypeI, vectorTypeS>::innerSolve(Eigen::VectorXd &rhs,
                                                     Eigen::VectorXd &result) {
  std::cout << "++++++++++++++++++ MKL: Solve *********************"
            << std::endl;
  /* -------------------------------------------------------------------- */
  /* .. Back substitution and iterative refinement. */
  /* -------------------------------------------------------------------- */
  double *x = (double *)mkl_calloc(rhs.size() * nrhs, sizeof(double), 64);
  phase = 33;
  iparm[7] = 0; /* Max numbers of iterative refinement steps. */

  PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &N_MKL, Ax, Ap, Ai, NULL, &nrhs,
          iparm, &msglvl, rhs.data(), x, &error);

  if (error != 0) {
    std::cerr << "ERROR during solve - code: " << error << std::endl;
  }
  result = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x, N_MKL);
  assert(Base::N == N_MKL);
  mkl_free(x);
  //    std::cout << "Res is: " << (rhs - A_CSC.selfadjointView<Eigen::Lower>()
  //    * result).norm() << std::endl;
}

// template <typename vectorTypeI, typename vectorTypeS>
// void MKLSolver<vectorTypeI, vectorTypeS>::updateMKL()
//{
//     Base::load_time = omp_get_wtime();
//
//     //    //    std::string f_name = "/home/behrooz/Desktop/oops.mtx";
//     //    //    std::ofstream f1_out(f_name);
//     //    //
//     //    A_CSC.conservativeResize(Base::numRows, Base::numRows);
//     //    A_CSC.reserve(Base::a.size());
//     //    //
//     //    std::copy(Base::ia.data(), Base::ia.data() + Base::ia.size(),
//     A_CSC.outerIndexPtr()); // pointer array (p)
//     //    std::copy(Base::ja.data(), Base::ja.data() + Base::ja.size(),
//     A_CSC.innerIndexPtr()); //  index array (i)
//     //    std::copy(Base::a.data(), Base::a.data() + Base::a.size(),
//     A_CSC.valuePtr()); // value array (x)
//
//     Ap = Base::ia.data();
//     Ai = Base::ja.data();
//     Ax = Base::a.data();
//     N = Base::numRows;
//     //    print_csc(N,
//     //        Ap,
//     //        Ai,
//     //        Ax,
//     //        f1_out.rdbuf(), "");
//     //    f1_out.close();
//
//     Base::load_time = omp_get_wtime() - Base::load_time;
//     Base::total_load_time += Base::load_time;
// }

template class MKLSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC

#endif