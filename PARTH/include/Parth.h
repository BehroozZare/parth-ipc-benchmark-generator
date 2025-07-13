//
// Created by behrooz zare on 2024-04-07.
//

#ifndef PARTH_SOLVER_PARTH_H
#define PARTH_SOLVER_PARTH_H

#include "HMD.h"
#include "Integrator.h"
#include "ParthTypes.h"
#include <tuple>
#include <vector>

namespace PARTH {
class Parth {
public:
  //============================ OPTIONS =================================
  ReorderingType reorder_type;
  int num_cores = 10;
  bool verbose = true;
  bool activate_aggressive_reuse = true;
  bool lagging = true;
  int ND_levels = 0;
  int sim_dim = 1;

  ///--------------------------------------------------------------------------
  /// setReorderingType - set the reordering type from METIS, AMD, AUTO
  ///--------------------------------------------------------------------------
  void setReorderingType(ReorderingType type);
  ReorderingType getReorderingType() const;

  ///--------------------------------------------------------------------------
  /// setVerbose - if true, report will be printed
  ///--------------------------------------------------------------------------
  void setVerbose(bool verbose);
  bool getVerbose() const;

  ///--------------------------------------------------------------------------
  /// setNumSubMesh - Set the number of nested dissection levels for HMD
  ///--------------------------------------------------------------------------
  void setNDLevels(int num);
  int getNDLevels() const;

  ///--------------------------------------------------------------------------
  /// setNumSubMesh - Set the number of nested dissection levels for HMD
  ///--------------------------------------------------------------------------
  void setRelocatedLevel(int lvl);
  int getRelocatedLevel() const;

  void setLagLevel(int lvl);
  int getLagLevel() const;


  ///--------------------------------------------------------------------------
  /// setNumberOfCores - Set the number of cores for region parallelism
  ///--------------------------------------------------------------------------
  void setNumberOfCores(int num_cores);
  int getNumberOfCores() const;


  //============================ VARIABLES ================================
  HMD hmd;
  Integrator integrator;

  int *Mp;
  int *Mi;
  int M_n;
  int M_n_prev;
  std::vector<int> Mp_prev;
  std::vector<int> Mi_prev;

  std::vector<int> mesh_perm;
  std::vector<std::pair<int, int>> changed_dof_edges;

  std::vector<int> new_to_old_map; // Map the id of new mesh dofs to old ones.
                                   // If the DOFs are new, the value is -1
  std::vector<int> old_to_new_map; // Map the id of old mesh dofs to new ones.
                                   // If the DOFs are removed, the value is -1
  std::vector<int> added_dofs;
  std::vector<int> removed_dofs;


  double dof_change_integrator_time = 0.0;
  double change_computation_time = 0.0;
  double map_mesh_to_matrix_computation_time = 0.0;
  double compute_permutation_time = 0.0;
  double init_time = 0.0;

  //============================ FUNCTIONS =================================
  Parth();
  ~Parth();

  void clearParth();

  ///--------------------------------------------------------------------------
  /// setMatrixPointers - it only set the pointers to the appropriate used
  /// defined allocated matrices\n NOTE: The correctness of memory allocation is
  /// on the user and Parth may change the arrays
  ///--------------------------------------------------------------------------
  void
  setMeshPointers(int n,                /// [in] Number of rows/columns
                  int *Mp,              ///<[in] pointer array
                  int *Mi,              ///<[in] index array
                  std::vector<int> &map ///<[in] the new to old dof index map
  );

  void setMeshPointers(int n,   /// [in] Number of rows/columns
                       int *Mp, ///<[in] pointer array
                       int *Mi  ///<[in] index array
  );

  ///--------------------------------------------------------------------------
  /// setDOFsMapper - Set the DOFs mapper where for each dof in the current mesh
  /// the corresponding if of the dof in previous mesh is stored. If the dof is
  /// new, the value is -1 -> Use this if you want reuse with adding and
  /// deleting DOFs
  ///--------------------------------------------------------------------------
  void setNewToOldDOFMap(std::vector<int> &map);
  void computeOldToNewDOFMap();

  /// computePermutation - Compute the permutation of the mesh using Parth
  /// algorithm
  void computePermutation(std::vector<int> &perm, int dim = 3);

  /// mapMeshPermToMatrixPerm - Map the mesh permutation to matrix permutation
  void mapMeshPermToMatrixPerm(std::vector<int> &mesh_perm,
                               std::vector<int> &matrix_perm, int dim = 3);

  double getReuse();
  int getNumChanges();
  void printTiming();
  void resetTimers();
};
} // namespace PARTH

#endif // PARTH_SOLVER_PARTH_H
