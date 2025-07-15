# Parth’s IPC Benchmark: Matrix Generation for Parth Evaluation

This repository is a fork of [IPC (Incremental Potential Contact)](https://github.com/ipc-sim/IPC), specifically modified for matrix generation. The objective of this fork is to generate sparse semi-positive definite matrices used to evaluate [Parth](https://arxiv.org/abs/2501.04011), which reuses symbolic analysis across sequences of calls to sparse solvers with dynamic sparsity patterns.

---

## Setup

### Installing IPC

To install IPC, please refer to the instructions in the [original IPC repository](https://github.com/ipc-sim/IPC). Follow their setup guide for your specific platform.

Howver, you can build this code-base using the standard CMake build process:

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j<#threads>
```

---

### Basic Benchmark Generation

After building the project, you can test whether your installation works correctly by running:

```bash
cd build
./IPC_bin 10 ../input/paperExamples/12_matOnBoard.txt
```

You should see the simulation of *matOnBoard* in the IPC benchmark. An output folder will be automatically created, containing Hessian matrices in the format:

```bash
hessian_<#frame>_<#iteration>_last_IPC.mtx
```

You can also try other examples in the `input/paperExamples` folder.

---

### Checkpoint and Resume Functionality

Simulations can be long-running, so IPC provides checkpoint functionality to resume from intermediate states. This is useful for:

1. **Parallelizing computation** across multiple machines
2. **Extracting matrices** from specific time points without running the full simulation
3. **Recovering from interruptions** without losing progress

The checkpoint system saves simulation state files that you can load to resume from a particular frame or iteration. This capability is leveraged to speed up benchmark generation.

---

### Fast Generation Scripts

In [EvaluationPipeline](EvaluationPipeline), you’ll find Python scripts that automate:

* Batch processing of multiple simulation segments
* Parallel execution across multiple cores or machines
* Matrix extraction at specific simulation steps

All evaluation scripts use our Apptainer container ([parth\_docker.def](EvaluationPipeline/parth_docker.def)), which includes all necessary libraries (such as Intel MKL). For macOS, you will need to run the code natively without the container, but the script generation can be easily adapted.

---

#### Informal Explanation

In practice, generating this benchmark from scratch takes considerable time. To make this process easier, I precomputed many simulations on a server over several days to gather a collection of checkpoint files. These files let you resume simulations without starting from the beginning. Specifically, after running each simulation end-to-end, I divided it into 10 segments that can be used for faster iterative development.

If you have enough storage (approximately 10TB in my case), you can also store all the Hessian matrices for downstream evaluation of different linear solvers using Parth’s codebase. The [EvaluationPipeline/README.md](EvaluationPipeline/README.md) includes detailed instructions on downloading these status files (covering every 10% of the simulation) and configuring the evaluation pipeline, along with examples of configuration files and job submission workflows.

---

### Mac Users: OpenMP Configuration

Mac users may need to manually link the OpenMP library in CMake. You can find the relevant configuration in `CMakeLists.txt` around lines 295–305:

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

Adjust the path `/usr/local/Cellar/libomp/18.1.4/` as needed to match your local OpenMP installation.
(*Note: If someone provides a cleaner way to set this via a variable, contributions are welcome!*)

---

## Matrix Output Format

The generated Hessian matrices are saved in Matrix Market format (`.mtx`), which is compatible with most sparse matrix libraries and analysis tools. The matrices represent the system’s stiffness matrix at each simulation step, providing the sparse semi-positive definite matrices needed for Parth evaluation.

---

## Parth Integration

This fork includes Parth integration for advanced sparse matrix analysis and reuse evaluation. The components in the `PARTH/` directory provide:

* Symbolic analysis reuse
* Dynamic sparsity pattern handling
* Performance evaluation metrics
* Matrix permutation and reordering analysis

---

## Output Structure

Simulation outputs include:

* **Hessian matrices**: Sparse matrices in `.mtx` format
* **Performance metrics**: CSV files containing timing and reuse statistics
* **Simulation state**: Checkpoint files for resuming simulations
* **Visualization data**: PNG images and GIF animations (if enabled)

---

## Citation

If you use this fork in your research, please cite both the original IPC paper and the Parth paper:

* IPC: [Incremental Potential Contact: Intersection- and Inversion-free, Large-Deformation Dynamics](https://doi.org/10.1145/3386569.3392412)
* Parth: [Parth: Reusing Symbolic Analysis in Sparse Solvers with Dynamic Sparsity Patterns](https://arxiv.org/abs/2501.04011)