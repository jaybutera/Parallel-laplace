import sys
import os

f_sec1 = \
'''#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM
#SBATCH --ntasks-per-node '''
f_sec2 = \
'''
#SBATCH -t 08:00:00
export OMP_NUM_THREADS='''

f_mpi_sec1 = \
'''#!/bin/bash
#SBATCH -N '''
f_mpi_sec2 = \
'''
#SBATCH -p RM
#SBATCH --ntasks-per-node '''

def njob(num_cores):
    file_lines = f_sec1 + str(num_cores) + f_sec2 + str(num_cores) + '\nmpirun -np 1 ./' + sys.argv[1]

    with open('lap_job.sh', 'w') as f:
        f.write(file_lines)

def job_nodes(num_nodes, num_cores):
    file_lines = f_mpi_sec1 + str(num_nodes) + f_mpi_sec2 + str(num_cores) + \
    f_sec2 + str(num_cores) + '\nmpirun -np ' + str(num_nodes) + ' ./' + sys.argv[1]

    with open('lap_job.sh', 'w') as f:
        f.write(file_lines)

def run_ntasks(n_times):
    for _ in range(n_times):
        os.system('sbatch lap_job.sh')

# Takes a list of all node numbers to run n jobs on
def full_mpi_test(nodes, n):
    for i in nodes:
        '''
        dirname = './slurm1000_n' + str(i)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        '''

        job_nodes(i, 1)
        run_ntasks(n)


#njob(4)
#run_ntasks(5)
if len(sys.argv) < 2:
    print 'Specify a program name dumbass.'
else:
    #job_nodes(3,1)
    njob(8)
    run_ntasks(5)
    #full_mpi_test([1,2,3,4,8,14],5)

'''
    argv[1] is program name
'''
