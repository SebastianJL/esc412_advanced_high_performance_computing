#!/bin/bash -l
#SBATCH --job-name="main_cuda"
#SBATCH --account="uzg2"
#SBATCH --mail-type=None
#SBATCH --mail-user=johannes.lade@uzh.ch
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=2
#SBATCH --partition=debug
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread
#SBATCH -o slurm/%j-mpi-openmp.out
#SBATCH -e slurm/%j-mpi-openmp.err


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun make-build/mpi_openmp /store/uzh/uzh8/ESC412/ic_512.std 512
