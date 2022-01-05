#!/bin/bash -l
#SBATCH --job-name="ex13-2"
#SBATCH --account="uzg2"
#SBATCH --mail-type=None
#SBATCH --mail-user=johannes.lade@uzh.ch
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread
#SBATCH -o slurm/ex13-1-out-%j


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

srun make-build/ex13-2 /store/uzh/uzh8/ESC412/ic_512.std 512

