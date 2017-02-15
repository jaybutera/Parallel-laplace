import sys
import os

f_sec1 = \
'''#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node '''
f_sec2 = \
'''
#SBATCH -t 00:05:00
export OMP_NUM_THREADS='''

def njob(num_cores):
	file_lines = f_sec1 + str(num_cores) + f_sec2 + str(num_cores) + '\n./' + sys.argv[1]

	with open('lap_job.sh', 'w') as f:
		f.write(file_lines)

def run_ntasks(n_times):
	for _ in range(n_times):
		os.system('sbatch lap_job.sh')

#njob(4)
run_ntasks(2)
