#!/bin/bash
#SBATCH -J apogeebh          # job name
#SBATCH -o apogeebh.o%j             # output file name (%j expands to jobID)
#SBATCH -e apogeebh.e%j             # error file name (%j expands to jobID)
#SBATCH -n 1
#SBATCH -t 06:00:00             # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=adrn@princeton.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

cd /tigress/adrianp/projects/apogeebh/scripts/

source activate twoface

date

python load_dr15_db.py --db=../cache/apogeebh.sqlite \ --allstar=../apogeebh/db/tests/allStar.fits \ --allvisit=../apogeebh/db/tests/allVisit.fits -o

date
