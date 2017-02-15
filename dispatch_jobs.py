f_sec1 = \
'''#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node '''
f_sec2 = \
'''
#SBATCH -t 00:05:00
export OMP_NUM_THREADS='''
f_sec3 = \
'''
./stripe
'''

def njob(num_cores):
	file_lines = f_sec1 + str(num_cores) + f_sec2 + str(num_cores) + f_sec3

	with open('lap_job.sh', 'w') as f:
		f.write(file_lines)

njob(4)
