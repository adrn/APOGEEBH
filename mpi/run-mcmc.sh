#!/bin/bash
#SBATCH -J mcmc          # job name
#SBATCH -o mcmc.o%j             # output file name (%j expands to jobID)
#SBATCH -e apogeebh.e%j             # error file name (%j expands to jobID)
#SBATCH -n 224
#SBATCH -t 02:00:00             # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=adrn@princeton.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

cd /tigress/adrianp/projects/apogeebh/scripts/

module load openmpi/gcc/1.10.2/64

source activate twoface

date

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M05335134+2625329 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M06185816+4002167 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M06403354+3441198 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M01231070+1801407 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M11424549+0052349 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M12474509+1230192 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M00444105+8351358 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M05121632+4558157 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M08195523+2859517 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M09304744+5417137 --mpi

srun python run_mcmc.py -vv -c ../config/bh.yml --data-path=../data/candidates/ --apogeeid=2M17212080+6003296 --mpi

date
