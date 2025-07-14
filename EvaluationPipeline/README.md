# EvaluationPipeline

This folder contains scripts for generating batch job scripts and configuration files for running large-scale simulation tests on HPC clusters (e.g., Niagara, Cedar). The scripts automate the creation of SLURM batch scripts and simulation configuration files for different solvers, parameter sweeps, and segmented runs.

## Scripts Overview

### createFullSimulationScripts.py

- **Purpose:**
  - Generates SLURM batch scripts for running full (non-segmented) simulations for a set of test cases and parameter sweeps.
  - Organizes output scripts and CSV folders for each simulation.
- **Usage:**
  - Edit the `simulations` dictionary and the main block to specify which simulations and parameters to run.
  - Run the script to generate a directory of batch scripts and a master `run_all.sh` script for submitting all jobs.

### createSegSimulationScript.py

- **Purpose:**
  - Generates SLURM batch scripts for running simulations in segments (splitting the total frames into multiple jobs).
  - Also generates segmented configuration files for each simulation segment.
- **Usage:**
  - Edit the `simulations` dictionary and the main block to specify simulations, segmentation, and parameters.
  - Run the script to generate segmented config files and batch scripts for each segment, plus a master `run_all.sh`.

### createCheckPointScript.py

- **Purpose:**
  - Similar to the other scripts, but focused on generating scripts for checkpoint/restart testing and analysis.
  - Supports different solvers and parameter sweeps for checkpointed runs.
- **Usage:**
  - Edit the `simulations` dictionary and the main block as needed.
  - Run the script to generate checkpoint test scripts and a master submission script.

### parth_docker.def

- **Purpose:**
  - Singularity definition file for building the container used in the simulation jobs.
- **Usage:**
  - Use with Singularity to build the required container image for the simulation environment.

## General Workflow

1. Edit the relevant script to specify your simulations, parameters, and cluster settings.
2. Run the script to generate the batch scripts and config files.
3. Submit the generated `run_all.sh` script on your cluster to queue all jobs.

---

**Author:** Behrooz Zarebavani (2025)
