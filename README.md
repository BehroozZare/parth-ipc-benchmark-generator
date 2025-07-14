# Parth’s IPC Benchmark: Matrix Generation for Parth Evaluation

This repository is a fork of [IPC (Incremental Potential Contact)](https://github.com/ipc-sim/IPC), specifically modified for matrix generation. The objective of this fork is to generate sparse semi-positive definite matrices used to evaluate [Parth](https://arxiv.org/abs/2501.04011), which reuses symbolic analysis across sequences of calls to sparse solvers with dynamic sparsity patterns.

## Setup

### Installing IPC

To install IPC, please refer to the instructions in the [original IPC repository](https://github.com/ipc-sim/IPC). Follow their setup guide for your specific platform.

### Docker/Apptainer Usage

For containerized execution, we provide Apptainer (formerly Singularity) containers. The container files are referenced in the evaluation scripts and can be used to ensure consistent execution environments across different systems. For an overview of the evaluation pipeline and the scripts used to generate batch jobs and configuration files, see [EvaluationPipeline/README.md](EvaluationPipeline/README.md).

### Mac Users: OpenMP Configuration

Mac users may need to manually link the OpenMP library in CMake. The relevant configuration is located in `CMakeLists.txt` around lines 295–305:

```cmake
if (APPLE)
    set(OMP_BREW /usr/local/Cellar/libomp/18.1.4/)
    target_include_directories(${PROJECT_NAME}_dev PUBLIC ${OMP_BREW}/include)
    find_library(OMP NAMES omp HINTS ${OMP_BREW}/lib)
    if (OMP)
        message(STATUS "Found omp: ${OMP}")
    else ()
        message(FATAL_ERROR "omp library not found")
    endif ()
    target_link_libraries(${PROJECT_NAME}_dev PUBLIC ${OMP})
    add_definitions(-DOPENMP)
endif ()
```

You may need to adjust the path `/usr/local/Cellar/libomp/18.1.4/` to match your OpenMP installation location. (I am very lazy, so if someone kindly provides a default variable to set, I would appreciate it.)

## Creating the Benchmark

To create the benchmark, run the scripts with CHOLMOD as the default solver. This saves the Hessians in the same folder where IPC stores status files and related data.

### Running Simulations and Saving Hessians

To run simulations and save the Hessian matrices, use the following command format:

```bash
./build/IPC_bin <mode> <input_file> <parameters> <threads>
```

For example:

```bash
./build/IPC_bin 100 input/paperExamples/12_matOnBoard.txt 0.999 666 4 t12
```

The Hessian matrices are automatically saved in Matrix Market format (`.mtx`) in the output directory using the naming convention:

```
hessian_<frame>_<iteration>_<tag>_IPC.mtx
```

### Checkpoint and Resume Functionality

Because simulations can be long-running, checkpoint functionality allows you to resume from intermediate states. This is useful for:

1. **Parallelizing simulation computation** across multiple machines
2. **Extracting matrices** from specific time points without running the full simulation
3. **Recovering from interruptions** without losing progress

The checkpoint system saves simulation state files that can be used to resume from specific frames or iterations.

### Generation Scripts

Several Jupyter notebooks in the `EvaluationScriptGenerator/` directory help you create simulation scripts:

* `createCheckPointScript.ipynb`: Generates scripts for checkpoint-based simulations
* `createFullSimulationScripts.ipynb`: Creates scripts for full simulations
* `createSegSimulationScript.ipynb`: Generates scripts for segmented simulations

These scripts can be customized for your server environment and simulation requirements. They handle:

* Batch processing of multiple input files
* Parallel execution across multiple cores or machines
* Automatic checkpoint creation and resume functionality
* Matrix extraction at specific simulation steps

### Batch Processing

To run multiple simulations, use the provided batch script:

```bash
python batch.py --offline
```

This processes all input files in the `input/` directory with the appropriate thread configurations.

## Matrix Output Format

The generated Hessian matrices are saved in Matrix Market format (`.mtx`), which is compatible with most sparse matrix libraries and analysis tools. The matrices represent the system’s stiffness matrix at each simulation step, providing the sparse semi-positive definite matrices needed for Parth evaluation.

## PARTH Integration

This fork includes PARTH integration for advanced sparse matrix analysis and reuse evaluation. The PARTH components, located in the `PARTH/` directory, provide functionality for:

* Symbolic analysis reuse
* Dynamic sparsity pattern handling
* Performance evaluation metrics
* Matrix permutation and reordering analysis

## Output Structure

Simulation outputs include:

* **Hessian matrices**: Sparse matrices in `.mtx` format
* **Performance metrics**: CSV files containing timing and reuse statistics
* **Simulation state**: Checkpoint files for resuming simulations
* **Visualization data**: PNG images and GIF animations (if enabled)

## Citation

If you use this fork in your research, please cite both the original IPC paper and the Parth paper:

* IPC: [Incremental Potential Contact: Intersection- and Inversion-free, Large-Deformation Dynamics](https://doi.org/10.1145/3386569.3392412)
* Parth: [Parth: Reusing Symbolic Analysis in Sparse Solvers with Dynamic Sparsity Patterns](https://arxiv.org/abs/2501.04011)