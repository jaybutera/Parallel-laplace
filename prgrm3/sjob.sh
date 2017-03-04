#!/bin/bash
#SBATCH -N 3
#SBATCH -p RM
#SBATCH --ntasks-per-node 1
#SBATCH -t 00:05:00
mpirun -np 3 ./a.out
