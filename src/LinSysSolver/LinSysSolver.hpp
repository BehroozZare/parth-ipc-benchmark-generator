//
//  LinSysSolver.hpp
//  IPC
//
//  Created by Minchen Li on 6/30/18.
//
#pragma once

#include "Types.hpp"
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Parth.h>
#include <csv_utils.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <unordered_map>
#ifdef USE_TBB
#include <tbb/parallel_for.h>
#endif

#define CONFIG_PAR_PARTH 0
#define CONFIG_PAR_ORDER 1
#define CONFIG_PAR_RES 2

namespace IPC {

enum class LinSysSolverType {
  CHOLMOD,
  AMGCL,
  EIGEN,
  MKL,
  STRUMPACK,
  CG
};

enum class FillReducingType {
  METIS,
  AMD,
};

template <typename vectorTypeI, typename vectorTypeS> class LinSysSolver {
protected:
  int numRows;
  Eigen::VectorXi ia,
      ja; // ia is the outer array and ja is the inner array in CSC format

  ///@brief A full CSC matrix with no diagonal entries
  int M_n;
  std::vector<int> Mp;
  std::vector<int> Mi;
  std::vector<int> elements_group;

  std::vector<std::map<int, int>> IJ2aI;
  Eigen::VectorXd a;

  FillReducingType fill_reducing_type = FillReducingType::METIS;

  double numerical_reuse = 0;
  double max_possible_numerical_reuse = 0;
  int dim = 3;

  std::string csv_address;
  bool active_parth;

public:
  PARTH::Parth parth;
  std::vector<int> perm;

  virtual ~LinSysSolver(void){};

  static LinSysSolver<vectorTypeI, vectorTypeS> *
  create(const LinSysSolverType type, const std::vector<double> &comp_par,
         const std::string &Comp = "NONE");
  virtual LinSysSolverType type() const = 0;

public:
  int frame;
  int iteration;
  double load_time = 0;
  double analyze_time = 0;
  double factor_time = 0;
  double solve_time = 0;
  double parth_time = 0;
  double residual_comp_time = 0;

  double lazy_analysis = 0;
  double lazy_factorization = 0;

  double total_load_time = 0;
  double total_factor_time = 0;
  double total_solve_time = 0;
  double total_analyze_time = 0;

  int L_NNZ = 0;
  int NNZ = 0;
  int N = 0;
  double residual = 0;
  double req_res = -1;
  std::string outputFolderPath;
  std::string csv_data_address;

  std::vector<int> &getMeshP() { return Mp; }

  std::vector<int> &getMeshI() { return Mi; }

  std::vector<int> &getElemGroups() { return elements_group; }

  void setElemGroups(std::vector<int> &groups) { elements_group = groups; }

  virtual void copyVNeighbor(const std::vector<std::set<int>> &vNeighbor,
                             const std::set<int> &fixedVert) {
    M_n = vNeighbor.size();
    Mp.clear();
    Mi.clear();
    Mp.resize(M_n + 1);
    int M_nnz = 0;
    for (int col = 0; col < vNeighbor.size(); col++) {
      Mp[col] = M_nnz;
      for (const auto &iter : vNeighbor[col]) {
        Mi.emplace_back(iter);
        M_nnz++;
      }
    }
    Mp[vNeighbor.size()] = M_nnz;
    assert(Mi.size() == M_nnz);
  }

  virtual void set_pattern(const std::vector<std::set<int>> &vNeighbor,
                           const std::set<int> &fixedVert) {
    copyVNeighbor(vNeighbor, fixedVert);
    numRows = static_cast<int>(vNeighbor.size()) * DIM;
    ia.resize(numRows + 1);
    ia[0] = 1;       // 1 + nnz above row i
    ja.resize(0);    // colI of each element
    IJ2aI.resize(0); // map from matrix index to ja index
    IJ2aI.resize(numRows);

    std::vector<Eigen::VectorXi> ja_v(vNeighbor.size());
    std::vector<int> rowNNZ(numRows);
#ifdef USE_TBB
    tbb::parallel_for(0, (int)vNeighbor.size(), 1,
                      [&](int vI) {
#else
    for (int vI = 0; vI < vNeighbor.size(); ++vI) {
#endif
                        ja_v[vI].resize((vNeighbor[vI].size() + 1) * DIM);

                        ja_v[vI][0] = vI * DIM;
                        ja_v[vI][1] = ja_v[vI][0] + 1;
                        IJ2aI[ja_v[vI][0]][ja_v[vI][0]] = 0;
                        IJ2aI[ja_v[vI][0]][ja_v[vI][1]] = 1;
                        if constexpr (DIM == 3) {
                          ja_v[vI][2] = ja_v[vI][0] + 2;
                          IJ2aI[ja_v[vI][0]][ja_v[vI][2]] = 2;
                        }

                        int nnz = DIM;
                        for (const auto &nbVI : vNeighbor[vI]) {
                          if (nbVI > vI) {
                            ja_v[vI][nnz] = nbVI * DIM;
                            ja_v[vI][nnz + 1] = ja_v[vI][nnz] + 1;
                            IJ2aI[ja_v[vI][0]][ja_v[vI][nnz]] = nnz;
                            IJ2aI[ja_v[vI][0]][ja_v[vI][nnz + 1]] = nnz + 1;
                            if constexpr (DIM == 3) {
                              ja_v[vI][nnz + 2] = ja_v[vI][nnz] + 2;
                              IJ2aI[ja_v[vI][0]][ja_v[vI][nnz + 2]] = nnz + 2;
                            }
                            nnz += DIM;
                          }
                        }

                        rowNNZ[ja_v[vI][0]] = nnz;
                        if constexpr (DIM == 2) {
                          ja_v[vI].conservativeResize(nnz * DIM - 1);
                          ja_v[vI].tail(nnz - 1) = ja_v[vI].segment(1, nnz - 1);

                          IJ2aI[ja_v[vI][0] + 1] = IJ2aI[ja_v[vI][0]];
                          IJ2aI[ja_v[vI][0] + 1].erase(ja_v[vI][0]);

                          rowNNZ[ja_v[vI][0] + 1] = nnz - 1;
                        } else {
                          ja_v[vI].conservativeResize(nnz * DIM - 3);
                          ja_v[vI].segment(nnz, nnz - 1) =
                              ja_v[vI].segment(1, nnz - 1);
                          ja_v[vI].tail(nnz - 2) = ja_v[vI].segment(2, nnz - 2);

                          IJ2aI[ja_v[vI][0] + 1] = IJ2aI[ja_v[vI][0]];
                          IJ2aI[ja_v[vI][0] + 1].erase(ja_v[vI][0]);
                          IJ2aI[ja_v[vI][0] + 2] = IJ2aI[ja_v[vI][0] + 1];
                          IJ2aI[ja_v[vI][0] + 2].erase(ja_v[vI][0] + 1);

                          rowNNZ[ja_v[vI][0] + 1] = nnz - 1;
                          rowNNZ[ja_v[vI][0] + 2] = nnz - 2;
                        }
                      }
#ifdef USE_TBB
    );
#endif

    for (int rowI = 0; rowI < numRows; ++rowI) {
      ia[rowI + 1] = ia[rowI] + rowNNZ[rowI];
    }

    ja.resize(ia[numRows] - 1);
#ifdef USE_TBB
    tbb::parallel_for(0, (int)vNeighbor.size(), 1,
                      [&](int vI) {
#else
    for (int vI = 0; vI < vNeighbor.size(); ++vI) {
#endif
                        int rowIStart = vI * DIM;

                        ja.segment(ia[rowIStart] - 1, ja_v[vI].size()) =
                            ja_v[vI];

                        for (auto &indexI : IJ2aI[rowIStart]) {
                          indexI.second += ia[rowIStart] - 1;
                        }
                        for (auto &indexI : IJ2aI[rowIStart + 1]) {
                          indexI.second += ia[rowIStart + 1] - 2;
                        }
                        if constexpr (DIM == 3) {
                          for (auto &indexI : IJ2aI[rowIStart + 2]) {
                            indexI.second += ia[rowIStart + 2] - 3;
                          }
                        }
                      }
#ifdef USE_TBB
    );
#endif
    ja.array() += 1;
    a.resize(ja.size());

    // NOTE: fixed verts nnz entries are not eliminated
  }

  virtual void load(const char *filePath, Eigen::VectorXd &rhs) {
    FILE *in = fopen(filePath, "rb");
    assert(in);

    size_t vecSize;
    fread(&vecSize, sizeof(size_t), 1, in);
    std::cout << "ia size " << vecSize << std::endl;
    ia.resize(vecSize);
    fread(ia.data(), sizeof(ia[0]), vecSize, in);

    fread(&vecSize, sizeof(size_t), 1, in);
    std::cout << "ja size " << vecSize << std::endl;
    ja.resize(vecSize);
    fread(ja.data(), sizeof(ja[0]), vecSize, in);

    if (ia[0] == 0) {
      ia.array() += 1;
      ja.array() += 1;
    }

    fread(&vecSize, sizeof(size_t), 1, in);
    std::cout << "a size " << vecSize << std::endl;
    a.resize(vecSize);
    fread(a.data(), sizeof(a[0]), vecSize, in);

    fread(&vecSize, sizeof(size_t), 1, in);
    std::cout << "rhs size " << vecSize << std::endl;
    rhs.resize(vecSize);
    fread(rhs.data(), sizeof(rhs[0]), vecSize, in);

    numRows = vecSize;

    fclose(in);
    std::cout << "load done" << std::endl;
  }
  virtual void write(const char *filePath, const Eigen::VectorXd &rhs) {
    FILE *out = fopen(filePath, "wb");

    size_t vecSize = ia.size();
    fwrite(&vecSize, sizeof(vecSize), 1, out);
    fwrite(ia.data(), sizeof(ia[0]), ia.size(), out);

    vecSize = ja.size();
    fwrite(&vecSize, sizeof(vecSize), 1, out);
    fwrite(ja.data(), sizeof(ja[0]), ja.size(), out);

    vecSize = a.size();
    fwrite(&vecSize, sizeof(vecSize), 1, out);
    fwrite(a.data(), sizeof(a[0]), a.size(), out);

    vecSize = rhs.size();
    fwrite(&vecSize, sizeof(vecSize), 1, out);
    fwrite(rhs.data(), sizeof(rhs[0]), rhs.size(), out);

    fclose(out);
  }

  virtual void set_pattern(const Eigen::SparseMatrix<double> &mtr) {
    // NOTE: mtr must be SPD

    numRows = static_cast<int>(mtr.rows());

    ja.conservativeResize(mtr.nonZeros());
    memcpy(ja.data(), mtr.innerIndexPtr(),
           mtr.nonZeros() * sizeof(mtr.innerIndexPtr()[0]));

    ia.conservativeResize(numRows + 1);
    memcpy(ia.data(), mtr.outerIndexPtr(),
           (numRows + 1) * sizeof(mtr.outerIndexPtr()[0]));

    a.conservativeResize(mtr.nonZeros());
    memcpy(a.data(), mtr.valuePtr(),
           mtr.nonZeros() * sizeof(mtr.valuePtr()[0]));
  }

  virtual void innerAnalyze_pattern() = 0;
  virtual void analyze_pattern(void) {
    analyze_time = omp_get_wtime();
    innerAnalyze_pattern();
    analyze_time = omp_get_wtime() - analyze_time;
  }

  virtual bool innerFactorize(void) = 0;
  virtual bool factorize(void) {
    factor_time = omp_get_wtime();
    bool succeeded = innerFactorize();
    factor_time = omp_get_wtime() - factor_time;
    if (!succeeded) {
      std::cout << "---- Factorization failed" << std::endl;
    }
    return succeeded;
  }

  virtual void innerSolve(Eigen::VectorXd &rhs, Eigen::VectorXd &result) = 0;
  virtual void solve(Eigen::VectorXd &rhs, Eigen::VectorXd &result) {
    residual_comp_time = 0;
    solve_time = omp_get_wtime();
    innerSolve(rhs, result);
    solve_time = omp_get_wtime() - solve_time;

    // Compute residual
    Eigen::SparseMatrix<double> mtr;
    mtr.resize(numRows, numRows);
    mtr.reserve(ja.size());
    memcpy(mtr.outerIndexPtr(), ia.data(), ia.size() * sizeof(ia[0]));
    memcpy(mtr.innerIndexPtr(), ja.data(), ja.size() * sizeof(ja[0]));
    memcpy(mtr.valuePtr(), a.data(), a.size() * sizeof(a[0]));
    residual = (rhs - mtr.selfadjointView<Eigen::Lower>() * result).norm();
    std::cout << "---- The linear solver residual is: " << residual
              << std::endl;
    addTotalTime();
  }

  virtual void setReorderingType(FillReducingType reorderingType) {
    this->fill_reducing_type = reorderingType;
    if (reorderingType == FillReducingType::METIS) {
      parth.setReorderingType(PARTH::ReorderingType::METIS);
      std::cout << "Using METIS" << std::endl;
    } else if (reorderingType == FillReducingType::AMD) {
      parth.setReorderingType(PARTH::ReorderingType::AMD);
      std::cout << "Using AMD" << std::endl;
    }
  }

  virtual void setParthActivation(bool activate) {
    active_parth = activate;
    if (active_parth) {
      std::cout << "Activating Parth" << std::endl;
    } else {
      std::cout << "Deactivating Parth" << std::endl;
    }
  }

  virtual void setSimulationDIM(int dim) {
    this->dim = dim;
    this->parth.sim_dim = dim;
  }

  virtual void multiply(const Eigen::VectorXd &x, Eigen::VectorXd &Ax) {
    assert(x.size() == numRows);
    assert(IJ2aI.size() == numRows);

    Ax.setZero(numRows);
    for (int rowI = 0; rowI < numRows; ++rowI) {
      for (const auto &colI : IJ2aI[rowI]) {
        Ax[rowI] += a[colI.second] * x[colI.first];
        if (rowI != colI.first) {
          Ax[colI.first] += a[colI.second] * x[rowI];
        }
      }
    }
  }

public:
  virtual void outputFactorization(const std::string &filePath) {
    assert(0 && "please implement!");
  }
  virtual double coeffMtr(int rowI, int colI) const {
    if (rowI > colI) {
      // return only upper right part for symmetric matrix
      int temp = rowI;
      rowI = colI;
      colI = temp;
    }
    assert(rowI < IJ2aI.size());
    const auto finder = IJ2aI[rowI].find(colI);
    if (finder != IJ2aI[rowI].end()) {
      return a[finder->second];
    } else {
      return 0.0;
    }
  }
  virtual void getCoeffMtr(Eigen::SparseMatrix<double> &mtr) const {
    mtr.resize(numRows, numRows);
    mtr.setZero();
    mtr.reserve(a.size() * 2 - numRows);
    for (int rowI = 0; rowI < numRows; rowI++) {
      for (const auto &colIter : IJ2aI[rowI]) {
        mtr.insert(rowI, colIter.first) = a[colIter.second];
        if (rowI != colIter.first) {
          mtr.insert(colIter.first, rowI) = a[colIter.second];
        }
      }
    }
  }
  virtual void getCoeffMtr_lower(Eigen::SparseMatrix<double> &mtr) const {
    assert(numRows > 0);

    mtr.conservativeResize(numRows, numRows);
    mtr.reserve(a.size());

    memcpy(mtr.innerIndexPtr(), ja.data(), ja.size() * sizeof(ja[0]));
    memcpy(mtr.outerIndexPtr(), ia.data(), ia.size() * sizeof(ia[0]));
    memcpy(mtr.valuePtr(), a.data(), a.size() * sizeof(a[0]));
  }
  virtual void getTriplets(const Eigen::VectorXi &nodeList,
                           std::vector<Eigen::Triplet<double>> &triplet) const {
    std::map<int, int> rowIMapper;
    for (int i = 0; i < nodeList.size(); ++i) {
      int startI = i * DIM;
      int startRowI = nodeList[i] * DIM;

      rowIMapper[startRowI] = startI;
      rowIMapper[startRowI + 1] = startI + 1;
      if constexpr (DIM == 3) {
        rowIMapper[startRowI + 2] = startI + 2;
      }
    }

    triplet.resize(0);
    for (int rowI = 0; rowI < numRows; rowI++) {
      auto rowIFinder = rowIMapper.find(rowI);
      for (const auto &colIter : IJ2aI[rowI]) {
        auto colIFinder = rowIMapper.find(colIter.first);
        if (rowIFinder != rowIMapper.end() && colIFinder != rowIMapper.end()) {
          triplet.emplace_back(rowIFinder->second, colIFinder->second,
                               a[colIter.second]);
          if (rowIFinder->second != colIFinder->second) {
            triplet.emplace_back(colIFinder->second, rowIFinder->second,
                                 a[colIter.second]);
          }
        }
      }
    }
  }

  virtual void setCoeff(int rowI, int colI, double val) {
    if (rowI <= colI) {
      assert(rowI < IJ2aI.size());
      const auto finder = IJ2aI[rowI].find(colI);
      assert(finder != IJ2aI[rowI].end());
      a[finder->second] = val;
    }
  }
  virtual void setCoeff(const LinSysSolver<vectorTypeI, vectorTypeS> *other,
                        double multiplier) {
    assert(numRows == other->numRows);
    assert(ja.size() == other->a.size());

    a = multiplier * other->a;
  }
  virtual void setZero(void) { a.setZero(); }
  virtual void setUnit_row(int rowI) {
    assert(numRows == IJ2aI.size());
    assert(rowI < numRows);
    for (const auto &colIter : IJ2aI[rowI]) {
      a[colIter.second] = (colIter.first == rowI);
    }
  }
  virtual void setUnit_row(int rowI, std::unordered_map<int, double> &rowVec) {
    assert(numRows == IJ2aI.size());
    assert(rowI < numRows);
    rowVec.clear();
    for (const auto &colIter : IJ2aI[rowI]) {
      rowVec[colIter.first] = a[colIter.second];
      a[colIter.second] = (colIter.first == rowI);
    }
  }
  virtual void setUnit_col(int colI, const std::set<int> &rowVIs) {
    assert(numRows == IJ2aI.size());
    assert(colI < numRows);
    for (const auto &rowVI : rowVIs) {
      for (int dimI = 0; dimI < DIM; ++dimI) {
        int rowI = rowVI * DIM + dimI;
        assert(rowI < numRows);
        if (rowI <= colI) {
          const auto finder = IJ2aI[rowI].find(colI);
          if (finder != IJ2aI[rowI].end()) {
            a[finder->second] = (rowI == colI);
          }
        }
      }
    }
  }
  virtual void setUnit_col_dim1(int colI, const std::set<int> &rowVIs) {
    assert(numRows == IJ2aI.size());
    assert(colI < numRows);
    for (const auto &rowI : rowVIs) {
      assert(rowI < numRows);
      if (rowI <= colI) {
        const auto finder = IJ2aI[rowI].find(colI);
        if (finder != IJ2aI[rowI].end()) {
          a[finder->second] = (rowI == colI);
        }
      }
    }
  }

  virtual void addCoeff(int rowI, int colI, double val) {
    if (rowI <= colI) {
      assert(rowI < IJ2aI.size());
      const auto finder = IJ2aI[rowI].find(colI);
      assert(finder != IJ2aI[rowI].end());
      a[finder->second] += val;
    }
  }
  virtual void precondition_diag(const Eigen::VectorXd &input,
                                 Eigen::VectorXd &output) {
    assert(numRows == input.size());
    output.resize(numRows);
    for (int rowI = 0; rowI < numRows; ++rowI) {
      const auto finder = IJ2aI[rowI].find(rowI);
      assert(finder != IJ2aI[rowI].end());
      output[rowI] = input[rowI] / a[finder->second];
    }
  }
  virtual void getMaxDiag(double &maxDiag) {
    maxDiag = -std::numeric_limits<double>::infinity();
    for (int rowI = 0; rowI < numRows; ++rowI) {
      const auto finder = IJ2aI[rowI].find(rowI);
      assert(finder != IJ2aI[rowI].end());
      if (maxDiag < a[finder->second]) {
        maxDiag = a[finder->second];
      }
    }
  }
  virtual void addCoeff(const LinSysSolver<vectorTypeI, vectorTypeS> *other,
                        double multiplier) {
    assert(numRows == other->numRows);
    if (a.size() == other->a.size()) {
      a += multiplier * other->a;
    } else {
      for (int rowI = 0; rowI < numRows; ++rowI) {
        for (const auto &colIter : other->IJ2aI[rowI]) {
          const auto finder = IJ2aI[rowI].find(colIter.first);
          if (finder != IJ2aI[rowI].end()) {
            a[finder->second] += multiplier * other->a[colIter.second];
          }
        }
      }
    }
  }

  virtual void saveHessian(std::string tag) {
    std::string uniq_f_name =
        std::to_string(frame) + "_" + std::to_string(iteration);
    std::string f_name =
        outputFolderPath + "/hessian_" + uniq_f_name + "_" + tag + "_IPC.mtx";

    std::ofstream f1_out(f_name);
    size_t n = getNumRows();
    int *Ap = get_ia().data();
    int *Ai = get_ja().data();
    double *Ax = get_a().data();
    std::streambuf *out = f1_out.rdbuf();
    const std::string indent = "";
    const std::string &beg = "%%MatrixMarket matrix coordinate real symmetric";

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

  virtual int getNumRows(void) const { return numRows; }
  virtual int getNumNonzeros(void) const { return a.size(); }

  virtual int getFactorNonzeros(void) { return 0; }

  virtual const std::vector<std::map<int, int>> &getIJ2aI(void) const {
    return IJ2aI;
  }
  virtual Eigen::VectorXi &get_ia(void) { return ia; }
  virtual Eigen::VectorXi &get_ja(void) { return ja; }
  virtual Eigen::VectorXd &get_a(void) { return a; }
  virtual const Eigen::VectorXd &get_a(void) const { return a; }

  virtual int get_numCol_1D(void) { return M_n; }
  virtual int *get_p_1D(void) { return Mp.data(); }
  virtual int *get_i_1D(void) { return Mi.data(); }

  double getResidual(void) {
    std::cout << "The residual is " << residual << std::endl;
    return residual;
  }

  double getNumericReusePercentage() { return numerical_reuse; }
  double getHessianInfNorm() {
    return *std::max_element(a.data(), a.data() + a.size());
  }
  double getAxInfNorm() { return 0; }

  // PARTH STATISTICS
  double getContactSize() { return parth.getNumChanges(); }
  double getReusePercentage() { return parth.getReuse(); }
  double getIntegrationTime() { return parth.dof_change_integrator_time; }
  double getChangeDetTime() { return parth.dof_change_integrator_time; }
  double getPermCompTime() { return parth.compute_permutation_time; }
  double getMapTime() { return parth.map_mesh_to_matrix_computation_time; }
  double getTotalTime() {
    return parth.init_time + parth.dof_change_integrator_time +
           parth.change_computation_time + parth.compute_permutation_time +
           parth.map_mesh_to_matrix_computation_time;
  }

  virtual void addCSVRecord(std::string csv_address, int frame, int iter,
                            double end_to_end_time) {
    std::vector<std::string> Runtime_headers;
    Runtime_headers.emplace_back("N");
    Runtime_headers.emplace_back("NNZ");
    Runtime_headers.emplace_back("L_NNZ");

    Runtime_headers.emplace_back("nthreads");

    // Iteration ids
    Runtime_headers.emplace_back("FrameNum");
    Runtime_headers.emplace_back("NewtonIter");

    // Solve Quality
    Runtime_headers.emplace_back("residual");

    // Timing Headers
    Runtime_headers.emplace_back("load_time");
    Runtime_headers.emplace_back("analyze_time");
    Runtime_headers.emplace_back("parth_time");
    Runtime_headers.emplace_back("factor_time");
    Runtime_headers.emplace_back("solve_time");

    // Parth Statistic
    Runtime_headers.emplace_back("contact_size");
    Runtime_headers.emplace_back("parth_reuse");
    Runtime_headers.emplace_back("barb_reuse");
    Runtime_headers.emplace_back("max_possible_reuse");
    Runtime_headers.emplace_back("parth_integration_time");
    Runtime_headers.emplace_back("parth_change_det_time");
    Runtime_headers.emplace_back("parth_perm_comp_time");
    Runtime_headers.emplace_back("parth_map_time");
    Runtime_headers.emplace_back("parth_total_time_time");
    Runtime_headers.emplace_back("parth_reordering");

    Runtime_headers.emplace_back("total_load_time");
    Runtime_headers.emplace_back("total_analyze_time");
    Runtime_headers.emplace_back("total_factor_time");
    Runtime_headers.emplace_back("total_solve_time");

    Runtime_headers.emplace_back("iter_time");

    std::string Data_name;
    profiling_utils::CSVManager runtime_csv(csv_address, "some address",
                                            Runtime_headers, false);

    runtime_csv.addElementToRecord(N, "N");
    runtime_csv.addElementToRecord(NNZ, "NNZ");
    runtime_csv.addElementToRecord(L_NNZ, "L_NNZ");
    runtime_csv.addElementToRecord(parth.getNumberOfCores(), "nthreads");

    // Iteration ids
    runtime_csv.addElementToRecord(frame, "FrameNum");
    runtime_csv.addElementToRecord(iter, "NewtonIter");

    runtime_csv.addElementToRecord(residual, "residual");

    // Timing Headers
    runtime_csv.addElementToRecord(load_time, "load_time");
    runtime_csv.addElementToRecord(analyze_time, "analyze_time");
    runtime_csv.addElementToRecord(parth_time, "parth_time");
    runtime_csv.addElementToRecord(factor_time, "factor_time");
    runtime_csv.addElementToRecord(solve_time, "solve_time");

    // 3 Region STATISICS
    runtime_csv.addElementToRecord(getContactSize(), "contact_size");
    runtime_csv.addElementToRecord(getReusePercentage(), "parth_reuse");
    runtime_csv.addElementToRecord(numerical_reuse, "barb_reuse");
    runtime_csv.addElementToRecord(max_possible_numerical_reuse,
                                   "max_possible_reuse");
    runtime_csv.addElementToRecord(getIntegrationTime(),
                                   "parth_integration_time");
    runtime_csv.addElementToRecord(getChangeDetTime(), "parth_change_det_time");
    runtime_csv.addElementToRecord(getPermCompTime(), "parth_perm_comp_time");
    runtime_csv.addElementToRecord(getMapTime(), "parth_map_time");
    runtime_csv.addElementToRecord(getTotalTime(), "parth_total_time_time");
    if (fill_reducing_type == FillReducingType::METIS)
      runtime_csv.addElementToRecord("metis", "parth_reordering");
    else if (fill_reducing_type == FillReducingType::AMD)
      runtime_csv.addElementToRecord("amd", "parth_reordering");

    runtime_csv.addElementToRecord(total_load_time, "total_load_time");
    runtime_csv.addElementToRecord(total_analyze_time, "total_analyze_time");
    runtime_csv.addElementToRecord(total_factor_time, "total_factor_time");
    runtime_csv.addElementToRecord(total_solve_time, "total_solve_time");
    runtime_csv.addElementToRecord(end_to_end_time, "iter_time");

    runtime_csv.addRecord();
    resetTimer();
    parth.resetTimers();
    numerical_reuse = 0;
    max_possible_numerical_reuse = 0;
    parth.integrator.reuse_ratio = 1;
    parth.changed_dof_edges.clear();
  }

  void resetTimer() {
    analyze_time = 0;
    factor_time = 0;
    solve_time = 0;
    load_time = 0;
    parth_time = 0;
  }

  void addTotalTime() {
    total_load_time += load_time;
    total_analyze_time += analyze_time;
    total_factor_time += factor_time;
    total_solve_time += solve_time;
  }

  void resetTotalTime() {
    total_load_time = 0;
    total_analyze_time = 0;
    total_factor_time = 0;
    total_solve_time = 0;
  }

  void initParth() {
    if (FillReducingType::METIS == fill_reducing_type) {
      parth.setReorderingType(PARTH::ReorderingType::METIS);
    } else if (FillReducingType::AMD == fill_reducing_type) {
      parth.setReorderingType(PARTH::ReorderingType::AMD);
    } else {
      std::cerr << "UNKNOWN ORDERING TYPE - Setting PARTH reordering to METIS"
                << std::endl;
      parth.setReorderingType(PARTH::ReorderingType::METIS);
    }

    parth.setVerbose(true);
    parth.setNDLevels(6);
    parth.setNumberOfCores(10);
    parth.sim_dim = DIM;
    parth.activate_aggressive_reuse = false;
    parth.lagging = false;
  }

  void initConfig(std::vector<double> &config) {
    if (config.size() > 1) {
      if (config[CONFIG_PAR_PARTH] == 1) {
        std::cout << "Using PARTH ...." << std::endl;
        active_parth = true;
      } else {
        std::cout << "Not using PARTH ...." << std::endl;
        active_parth = false;
      }
    }

    if (config.size() > 2) {
      if (config[CONFIG_PAR_ORDER] == 0) {
        std::cout << "Using METIS ...." << std::endl;
        fill_reducing_type = FillReducingType::METIS;
        parth.setReorderingType(PARTH::ReorderingType::METIS);
      } else if (config[CONFIG_PAR_ORDER] == 1) {
        std::cout << "Using AMD ...." << std::endl;
        fill_reducing_type = FillReducingType::AMD;
        parth.setReorderingType(PARTH::ReorderingType::AMD);
      } else {
        std::cerr << "UNKNOWN ORDERING TYPE - Setting PARTH reordering to METIS"
                  << std::endl;
      }
    }

    if (config.size() > 3) {
      if (config[CONFIG_PAR_RES] != -1) {
        std::cout << "Setting residual to " << config[CONFIG_PAR_RES]
                  << std::endl;
        req_res = config[CONFIG_PAR_RES];
      }
    }

    initParth();
  }
};

} // namespace IPC
