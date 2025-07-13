#include "LinSysSolver.hpp"
#include "AMGCLSolver.hpp"
#include "CGSolver.hpp"
#include "CHOLMODSolver.hpp"
#include "EigenLibSolver.hpp"

#ifdef IPC_WITH_STRUMPACK
#include "STRUMPACKSolver.hpp"
#endif

#ifdef IPC_WITH_MKL
#include "MKLSolver.hpp"
#endif

#include <spdlog/spdlog.h>

namespace IPC {

template <typename vectorTypeI, typename vectorTypeS>
LinSysSolver<vectorTypeI, vectorTypeS> *
LinSysSolver<vectorTypeI, vectorTypeS>::create(
    const LinSysSolverType type, const std::vector<double> &comp_par,
    const std::string &Comp) {

  switch (type) {
#ifdef IPC_WITH_CHOLMOD
  case LinSysSolverType::CHOLMOD:
    return new CHOLMODSolver<vectorTypeI, vectorTypeS>(comp_par);
#endif

#ifdef IPC_WITH_STRUMPACK
  case LinSysSolverType::STRUMPACK:
    return new STRUMPACKSolver<vectorTypeI, vectorTypeS>(Comp, comp_par);
#endif

#ifdef IPC_WITH_MKL
  case LinSysSolverType::MKL:
    return new MKLSolver<vectorTypeI, vectorTypeS>(comp_par);
#endif

#ifdef IPC_WITH_AMGCL
  case LinSysSolverType::AMGCL:
    return new AMGCLSolver<vectorTypeI, vectorTypeS>();
#endif

  case LinSysSolverType::EIGEN:
    return new EigenLibSolver<vectorTypeI, vectorTypeS>();
  case LinSysSolverType::CG:
    return new CGSolver<vectorTypeI, vectorTypeS>();
  default:
    spdlog::error("Uknown linear system solver type: {}", type);
    throw fmt::format("Uknown linear system solver type: {}", type);
  }
}

template class LinSysSolver<Eigen::VectorXi, Eigen::VectorXd>;

} // namespace IPC
