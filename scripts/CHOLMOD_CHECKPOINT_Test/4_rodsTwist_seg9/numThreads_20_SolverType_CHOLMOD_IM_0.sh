#!/bin/bash



#SBATCH --account=rrg-mmehride
#SBATCH --job-name="IPC_4_rodsTwist_seg9_numThreads_20_SolverType_CHOLMOD_IM_0"
#SBATCH --output="IPC_4_rodsTwist_seg9_numThreads_20_SolverType_CHOLMOD_IM_0_%j.out"
#SBATCH --error="IPC_4_rodsTwist_seg9_numThreads_20_SolverType_CHOLMOD_IM_0_%j.error"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --export=ALL
#SBATCH -t 24:00:00
#SBATCH --constraint=cascade



export SINGULARITY_BIND="$SCRATCH/test/parth-ipc-benchmark-generator/:/mnt"
export SING_SIF=/"$SCRATCH/test/parth-ipc-benchmark-generator/EvaluationPipeline/parth_docker.sif"
export num_threads=20
export MKL_NUM_THREADS=$num_threads
export OMP_NUM_THREADS=$num_threads
export VECLIB_MAXIMUM_THREADS=$num_threads
export PROG_PATH=/mnt/build/IPC_bin
export Config=/mnt/input/seg/4_rodsTwist_seg9.txt
export progMode=100


singularity exec $SING_SIF bash -c "$PROG_PATH $progMode $Config --SimName=../../csv/4_rodsTwist_seg9_numThreads_20_SolverType_CHOLMOD_IM_0 --output=/mnt/CHOLMOD_CHECKPOINT_Test/4_rodsTwist_seg9/numThreads_20_SolverType_CHOLMOD_IM_0 --DoAnalysis=0 --SolverType=CHOLMOD --numThreads=20 --IM=0"
