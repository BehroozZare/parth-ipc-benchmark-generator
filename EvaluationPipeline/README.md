# EvaluationPipeline

This folder contains scripts for generating batch job scripts and configuration files for running large-scale simulation tests on HPC clusters (e.g., Niagara, Cedar). The scripts automate the creation of SLURM batch scripts and simulation configuration files for different solvers, parameter sweeps, and segmented runs.

## Scripts Overview

### createFullSimulationScripts.py
- **Purpose:**
  - Generates SLURM batch scripts for running full (non-segmented) simulations for a set of test cases and parameter sweeps.
  - Organizes output scripts and CSV folders for each simulation.
- **Usage:**
  1. Edit the `simulations` dictionary and the main block to specify which simulations and parameters to run.
  2. Run the script:
     ```bash
     python createFullSimulationScripts.py
     ```
  3. Submit all jobs with the generated master script:
     ```bash
     cd scripts/<OUTPUT_DIR>
     bash run_all.sh
     ```

### createSegSimulationScript.py
- **Purpose:**
  - Generates SLURM batch scripts for running simulations in segments (splitting the total frames into multiple jobs).
  - Also generates segmented configuration files for each simulation segment.
  - Note: It requires the status files to resume the run.
- **Usage:**
  1. Edit the `simulations` dictionary and the main block to specify simulations, segmentation, and parameters.
  2. Run the script:
     ```bash
     python createSegSimulationScript.py
     ```
  3. Submit all jobs with the generated master script:
     ```bash
     cd scripts/<OUTPUT_DIR>
     bash run_all.sh
     ```

### createCheckPointScript.py
- **Purpose:**
  - Similar to the other scripts, but focused on generating scripts for checkpoint/restart testing and analysis.
  - Supports different solvers and parameter sweeps for checkpointed runs.
- **Usage:**
  1. Edit the `simulations` dictionary and the main block as needed.
  2. Run the script:
     ```bash
     python createCheckPointScript.py
     ```
  3. Submit all jobs with the generated master script:
     ```bash
     cd scripts/<OUTPUT_DIR>
     bash run_all.sh
     ```

### parth_docker.def
- **Purpose:**
  - Apptainer definition file for building the container used in the simulation jobs.
  - Use the following command to build your Docker in a Linux environment.
  - To build the container image, run:
    ```bash
    apptainer build parth_docker.sif parth_docker.def
    ```
- **Usage:**
  - Use with Apptainer/Singularity to build the required container image for the simulation environment.

## General Workflow
1. Edit the relevant script to specify your simulations, parameters, and cluster settings.
2. Run the script to generate the batch scripts and config files (see above for commands).
3. Submit the generated `run_all.sh` script on your cluster to queue all jobs:
   ```bash
   cd scripts/<OUTPUT_DIR>
   bash run_all.sh
   ```

---
**Author:** Behrooz Zarebavani (2025)
