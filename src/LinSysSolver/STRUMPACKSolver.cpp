//
//  CHOLMODSolver.cpp
//  IPC
//
//  Created by Minchen Li on 6/22/18.
//

#ifdef IPC_WITH_STRUMPACK

#include "STRUMPACKSolver.hpp"
#include "getRSS.hpp"

#include <StrumpackSparseSolver.hpp>
#include <iostream>
#include <stdexcept>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
STRUMPACKSolver<vectorTypeI, vectorTypeS>::STRUMPACKSolver(std::string Comp,
    std::vector<double> strum_par)
{

    std::cout << "The compression format is " << Comp << std::endl;
    // Options
    this->Compression = Comp;

    if (!strum_par.empty()) {

        if (strum_par.size() >= 1) {
            sp.options().set_rel_tol(strum_par[0]);
        }

        if (strum_par.size() >= 2) {
            sp.options().set_abs_tol(strum_par[1]);
        }

        if (strum_par.size() >= 3) {
            if (strum_par[2] != 0) {
                sp.options().HSS_options().set_leaf_size(64);
            }
        }

        if (strum_par.size() >= 4) {
            if (strum_par[3] != 0) {
                sp.options().HSS_options().set_clustering_algorithm(strumpack::ClusteringAlgorithm::KD_TREE);
            }
        }

        if (strum_par.size() >= 5) {
            if (strum_par[4] == 0) {
                std::cout << "=========================== CHOOSING METIS REORDERING ===========================" << std::endl;
                sp.options().set_reordering_method(strumpack::ReorderingStrategy::METIS);
            }
            else if (strum_par[4] == 1) {
                std::cout << "=========================== CHOOSING SCOTCH REORDERING ===========================" << std::endl;
                sp.options().set_reordering_method(strumpack::ReorderingStrategy::SCOTCH);
            }
            else if (strum_par[4] == 2) {
                std::cout << "=========================== CHOOSING NATURAL REORDERING ===========================" << std::endl;
                sp.options().set_reordering_method(strumpack::ReorderingStrategy::NATURAL);
            }
            else {
                std::cout << "=========================== Invalid Reordering ===========================" << std::endl;
            }
        }

        if (strum_par.size() >= 6) {
            double curr_size = sp.options().compression_min_sep_size();
            std::cout << "The prev compression sep size was: " << curr_size << std::endl;
            if (strum_par[5] != 0) {
                std::cout << "=========================== Compression Separation Size:" << strum_par[5] * curr_size << "  ===========================" << std::endl;
                sp.options().set_compression_min_sep_size(strum_par[5] * curr_size);
            }
        }

        if (strum_par.size() >= 7) {
            double curr_size = sp.options().compression_min_front_size();
            std::cout << "The prev compression front size was: " << curr_size << std::endl;
            if (strum_par[6] != 0) {
                std::cout << "=========================== Compression Front Size:" << strum_par[6] * curr_size << "  ===========================" << std::endl;
                sp.options().compression_min_front_size(strum_par[6] * curr_size);
            }
        }

        if (strum_par.size() >= 8) {
            if (strum_par[7] != 0) {
                if (Comp == "HSS") {
                    double curr_tol = sp.options().HSS_options().rel_tol();
                    std::cout << "=========================== Compression rel_tol for HSS: " << strum_par[7] * curr_tol << "  ===========================" << std::endl;
                    std::cout << "The prev compression rel_tol: " << curr_tol << std::endl;
                    sp.options().HSS_options().set_rel_tol(strum_par[7] * curr_tol);
                }
                else if (Comp == "BLR") {
                    double curr_tol = sp.options().BLR_options().rel_tol();
                    std::cout << "=========================== Compression rel_tol for BLR: " << strum_par[7] * curr_tol << "  ===========================" << std::endl;
                    std::cout << "The prev compression rel_tol: " << curr_tol << std::endl;
                    sp.options().BLR_options().set_rel_tol(strum_par[7] * curr_tol);
                }
                else if (Comp == "None") {
                    std::cout << "Compression is not going to happen in None compression" << std::endl;
                }
                else {
                    std::cerr << "Unknown compression " << Comp << std::endl;
                }
            }
        }

        if (strum_par.size() >= 9) {
            if (strum_par[8] != 0) {
                if (Comp == "HSS") {
                    double curr_tol = sp.options().HSS_options().abs_tol();
                    std::cout << "=========================== Compression abs_tol for HSS: " << strum_par[8] * curr_tol << "  ===========================" << std::endl;
                    std::cout << "The prev compression rel_tol: " << curr_tol << std::endl;
                    sp.options().HSS_options().set_abs_tol(strum_par[8] * curr_tol);
                }
                else if (Comp == "BLR") {
                    double curr_tol = sp.options().BLR_options().abs_tol();
                    std::cout << "=========================== Compression abs_tol for BLR: " << strum_par[8] * curr_tol << "  ===========================" << std::endl;
                    std::cout << "The prev compression abs_tol: " << curr_tol << std::endl;
                    sp.options().BLR_options().set_abs_tol(strum_par[8] * curr_tol);
                }
                else if (Comp == "None") {
                    std::cout << "Compression is not going to happen in None compression" << std::endl;
                }
                else {
                    std::cerr << "Unknown compression " << Comp << std::endl;
                }
            }
        }
    }

    sp.options().set_gmres_restart(30); // ...
    sp.options().set_matching(strumpack::MatchingJob::NONE);
    if (Comp == "HSS") {
        sp.options().set_compression(strumpack::CompressionType::HSS); // enable HSS compression
    }
    else if (Comp == "BLR") {
        sp.options().set_compression(strumpack::CompressionType::BLR); // enable BLR compression
    }
    else if (Comp == "NONE") {
        sp.options().set_compression(strumpack::CompressionType::NONE); // no compression - exact solve
    }
    else {
        std::cerr << "UNKNOWN Compression" << std::endl;
        throw std::runtime_error("ERROR");
    }

    sp.options().set_verbose(true);
}

template <typename vectorTypeI, typename vectorTypeS>
void STRUMPACKSolver<vectorTypeI, vectorTypeS>::set_pattern(const std::vector<std::set<int>>& vNeighbor,
    const std::set<int>& fixedVert)
{
    Base::set_pattern(vNeighbor, fixedVert);

    Base::ia.array() -= 1;
    Base::ja.array() -= 1; // CHOLMOD's index starts from 0
    updateSP();

    //    int* row_ptr = Base::ia.data(); // N+1 integers
    //    int* col_ind = Base::ja.data(); // nnz integers
    //    double* val = Base::a.data(); // nnz scalars

    //    sp.set_csr_matrix(Base::numRows, row_ptr, col_ind, val, true); // set the matrix (copy)
    //    A->i = Base::ja.data();
    //    A->p = Base::ia.data();
    //    A->x = Base::a.data();
}
template <typename vectorTypeI, typename vectorTypeS>
void STRUMPACKSolver<vectorTypeI, vectorTypeS>::set_pattern(const Eigen::SparseMatrix<double>& mtr)
{
    Base::set_pattern(mtr);
    Base::numRows = static_cast<int>(mtr.rows());
    updateSP();
}

template <typename vectorTypeI, typename vectorTypeS>
void STRUMPACKSolver<vectorTypeI, vectorTypeS>::analyze_pattern(void)
{
    std::cout << "++++++++++++++++++ Strumpack: Analysing *********************" << std::endl;
    Base::analyze_time = omp_get_wtime();
    sp.reorder(); // reorder matrix
    Base::analyze_time = omp_get_wtime() - Base::analyze_time;
    Base::total_analyze_time += Base::analyze_time;
}

template <typename vectorTypeI, typename vectorTypeS>
bool STRUMPACKSolver<vectorTypeI, vectorTypeS>::factorize(void)
{
    std::cout << "++++++++++++++++++ Strumpack: Factorizing *********************" << std::endl;
    updateSP();
    Base::factor_time = omp_get_wtime();
    sp.factor();
    Base::factor_time = omp_get_wtime() - Base::factor_time;
    Base::total_factor_time += Base::factor_time;
    return true; // TODO:CHECK FOR SPD FLAGS LATER
}

template <typename vectorTypeI, typename vectorTypeS>
void STRUMPACKSolver<vectorTypeI, vectorTypeS>::solve(Eigen::VectorXd& rhs,
    Eigen::VectorXd& result)
{
    std::cout << "++++++++++++++++++ Strumpack: Solve *********************" << std::endl;
    result.conservativeResize(rhs.size());
    Base::solve_time = omp_get_wtime();
    sp.solve(rhs.data(), result.data(), true);
    Base::solve_time = omp_get_wtime() - Base::solve_time;
    Base::total_solve_time += Base::solve_time;
    residual = (rhs - H_FULL_CSR * result).norm();
}

template <typename vectorTypeI, typename vectorTypeS>
void STRUMPACKSolver<vectorTypeI, vectorTypeS>::updateSP()
{
    Base::load_time = omp_get_wtime();

    Eigen::SparseMatrix<double> H0;
    H0.conservativeResize(Base::numRows, Base::numRows);
    H0.reserve(Base::a.size());

    std::copy(Base::ia.data(), Base::ia.data() + Base::ia.size(), H0.outerIndexPtr()); // pointer array (p)
    std::copy(Base::ja.data(), Base::ja.data() + Base::ja.size(), H0.innerIndexPtr()); //  index array (i)
    std::copy(Base::a.data(), Base::a.data() + Base::a.size(), H0.valuePtr()); // value array (x)
    H_FULL_CSR = H0.selfadjointView<Eigen::Lower>();
    int N = H_FULL_CSR.rows(); // construct an NxN CSR matrix with nnz nonzeros
    int* row_ptr = H_FULL_CSR.outerIndexPtr(); // N+1 integers
    int* col_ind = H_FULL_CSR.innerIndexPtr(); // nnz integers
    double* val = H_FULL_CSR.valuePtr(); // nnz scalars
    sp.set_csr_matrix(N, row_ptr, col_ind, val, true); // set the matrix (copy)

    Base::load_time = omp_get_wtime() - Base::load_time;
    Base::total_load_time += Base::load_time;
}

template <typename vectorTypeI, typename vectorTypeS>
int STRUMPACKSolver<vectorTypeI, vectorTypeS>::getFactorNonzeros(void)
{
    return sp.factor_nonzeros();
}

template <typename vectorTypeI, typename vectorTypeS>
double STRUMPACKSolver<vectorTypeI, vectorTypeS>::getResidual()
{
    std::cout << "The residual is " << residual << std::endl;
    return residual;
}

template class STRUMPACKSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC

#endif
