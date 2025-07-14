# EvaluationPipeline

This folder contains scripts for generating batch job scripts and configuration files for running large-scale IPC simulations for benchmark generation on HPC clusters (e.g., Niagara, Cedar). The scripts automate the creation of SLURM batch scripts and simulation configuration files for different solvers, parameter sweeps, and segmented runs.

## Scripts Overview

### createFullSimulationScripts.py

* **Purpose:**

  * Generates SLURM batch scripts for running full (non-segmented) simulations for a set of test cases and parameter sweeps.
  * Organizes output scripts and CSV folders for each simulation.
* **Usage:**

  1. Edit the `simulations` dictionary and the main block to specify which simulations and parameters to run.
  2. Run the script:

     ```bash
     python createFullSimulationScripts.py
     ```
  3. Submit all jobs using the generated master script:

     ```bash
     cd scripts/<OUTPUT_DIR>
     bash run_all.sh
     ```

### createSegSimulationScript.py

* **Purpose:**

  * Generates SLURM batch scripts for running simulations in segments (splitting the total frames into multiple jobs).
  * Also generates segmented configuration files for each simulation segment.
  * **Note:** This requires status files to resume the run.
* **Usage:**

  1. Edit the `simulations` dictionary and the main block to specify simulations, segmentation, and parameters.
  2. Run the script:

     ```bash
     python createSegSimulationScript.py
     ```
  3. Submit all jobs using the generated master script:

     ```bash
     cd scripts/<OUTPUT_DIR>
     bash run_all.sh
     ```

### createCheckPointScript.py

* **Purpose:**

  * Similar to the other scripts, but focused on generating scripts for checkpoint/restart testing and analysis.
  * Supports different solvers and parameter sweeps for checkpointed runs.
* **Usage:**

  1. Edit the `simulations` dictionary and the main block as needed.
  2. Run the script:

     ```bash
     python createCheckPointScript.py
     ```
  3. Submit all jobs using the generated master script:

     ```bash
     cd scripts/<OUTPUT_DIR>
     bash run_all.sh
     ```

### parth\_docker.def

* **Purpose:**

  * Apptainer definition file for building the container used in simulation jobs.
  * To build the container image, run:

    ```bash
    apptainer build parth_docker.sif parth_docker.def
    ```
* **Usage:**

  * Use with Apptainer/Singularity to build the required container image for the simulation environment.

### build\_ipc.sh

* **Purpose:**

  * Script to build the IPC simulation binary inside the Apptainer/Singularity container.
  * Installs dependencies and compiles the IPC codebase using CMake and Make.
* **Usage:**

  1. Place `build_ipc.sh` inside the container or bind-mount it into the container.
  2. Run the build script:

     ```bash
     bash build_ipc.sh
     ```
  3. The compiled binary (e.g., `IPC_bin`) will be placed in the specified build directory (e.g., `build/`).
* **Notes:**

  * Make sure the source code is available inside the container at the expected path.
  * You may need to adjust paths in the script to match your project structure.

## General Workflow

1. **If you are not running the full simulation from the beginning:**
   Download the checkpoints from [Zenodo](https://zenodo.org/records/15882563?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjFhMGZiMzM5LWQwN2ItNDZmMy1iNmZlLTFmOWM0NmI4ZGRkYiIsImRhdGEiOnt9LCJyYW5kb20iOiJkMzA3NTVmNGM4MTk0Zjg4ZGFiNzQyOTczM2MyYzI1ZiJ9.PZuqBIvVZMlDz-9r6Rwo28jFsSUTWxO3IDK3BjR7Ll2Y2zFb0oEdV6ouhIA7I3sGiJb_hg5LbJioRIUttWvXSw) or use your own checkpoint files and place them in the desired folder. (The default is a folder named `Checkpoint` In the root of this code base.)
   **Note:** If the checkpoints are not in the correct location, IPC will start the simulation from the beginning.
2. Edit the relevant script to specify your simulations, parameters, and cluster settings.
3. Run the script to generate the batch scripts and configuration files (see the usage sections above).
4. Submit the generated `run_all.sh` script on your cluster to queue all jobs:

   ```bash
   cd scripts/<OUTPUT_DIR>
   bash run_all.sh
   ```

---

**Author:** Behrooz Zarebavani (2025
