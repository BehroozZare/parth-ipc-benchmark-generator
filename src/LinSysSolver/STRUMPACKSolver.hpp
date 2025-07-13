//
//  MKL.hpp
//  IPC
//
//  Created by Behrooz
//
#pragma once

#ifndef IPC_STRUMPACK_SOLVER
#define IPC_STRUMPACK_SOLVER

#ifdef IPC_WITH_STRUMPACK

#include "LinSysSolver.hpp"
#include "omp.h"
#include "cholmod.h"

#include <Eigen/Eigen>

#include <vector>
#include <set>

#include <StrumpackSparseSolver.hpp>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
class STRUMPACKSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
    typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
    strumpack::StrumpackSparseSolver<double, int> sp; // create solver object
    Eigen::SparseMatrix<double> H_FULL_CSR;
    std::string Compression; // Options
    double residual;

public:
    STRUMPACKSolver(std::string Comp,
        std::vector<double> strum_par);
    ~STRUMPACKSolver(void) = default;

    LinSysSolverType type() const override { return LinSysSolverType::STRUMPACK; }

    void set_pattern(const std::vector<std::set<int>>& vNeighbor, const std::set<int>& fixedVert) override;
    void set_pattern(const Eigen::SparseMatrix<double>& mtr) override; // NOTE: mtr must be SPD

    void analyze_pattern(void) override;

    bool factorize(void) override;

    void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result) override;

    void updateSP();

    int getFactorNonzeros(void) override;
    double getResidual(void) override;

    //    virtual void multiply(const Eigen::VectorXd& x,
    //        Eigen::VectorXd& Ax) override;
    //
    //    virtual void outputFactorization(const std::string& filePath) override;
};

} // namespace IPC

#endif
#endif
