#!/bin/bash -l
#SBATCH --job-name="main_cuda"
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
#SBATCH -o slurm/cuda.out
#SBATCH -e slurm/cuda.err


export CRAY_CUDA_MPS=1

#srun make-build/cuda /store/uzh/uzh8/ESC412/ic_512.std 512
srun make-build/cuda ./input/b0-final.std 64
