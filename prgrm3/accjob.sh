#!/bin/bash
#SBATCH -N 3
#SBATCH -p GPU
#SBATCH -t 00:05:00
#SBATCH --gres=gpu:p100:2
mpirun -np 3 ./acc
