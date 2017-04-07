#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU-shared
#SBATCH --ntasks-per-node 1
#SBATCH --gres=gpu:p100:2
#SBATCH -t 00:10:00
./a.out
