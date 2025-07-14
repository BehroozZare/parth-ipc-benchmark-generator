#!/bin/bash
# A one-to-one conversion of your Singularity build script to Apptainer

set -euo pipefail            # good practice: fail fast & echo unset vars

#--- Bind mounts -------------------------------------------------------------
# Apptainer’s binding path
# export APPTAINER_BIND="Address/To/IPC/Cloned/Repository:/mnt"   # :contentReference[oaicite:0]{index=0}
export APPTAINER_BIND="$SCRATCH/test/parth-ipc-benchmark-generator:/mnt"


#--- Container image ---------------------------------------------------------
# Apptainer doesn’t care what you call this variable, but for clarity we
# rename it.  If you prefer the old name simply keep using $SING_SIF.
APP_SIF="${APP_SIF:-$SCRATCH/test/parth-ipc-benchmark-generator/EvaluationPipeline/parth_docker.sif}"

#--- Build inside the container ---------------------------------------------
apptainer exec "${APP_SIF}" bash -c '
    cd /mnt/                    && \
    mkdir -p build              && \
    cd build                    && \
    rm -f CMakeCache.txt        && \
    export STRUMPACKROOT=/mnt/STRUMPACK-6.3.1/install && \
    pwd                         && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j10 IPC_bin
'
