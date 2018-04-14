#!/bin/bash
#SBATCH -J apogeebh          # job name
#SBATCH -o apogeebh.o%j             # output file name (%j expands to jobID)
#SBATCH -e apogeebh.e%j             # error file name (%j expands to jobID)
#SBATCH -n 224
#SBATCH -t 12:00:00             # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=adrn@princeton.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

cd /tigress/adrianp/projects/apogeebh/scripts/

module load openmpi/gcc/1.10.2/64

source activate twoface

date

srun python run_bh.py -v -c ../config/bh.yml  --mpi --data-path=../data/candidates/ --ext=ecsv

date
