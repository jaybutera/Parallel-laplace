#!/bin/bash
#SBATCH -N 1
#SBATCH -p GPU
#SBATCH -t 00:05:00
#SBATCH --gres=gpu:p100:2
./acc
