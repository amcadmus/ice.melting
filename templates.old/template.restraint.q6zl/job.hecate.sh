#!/bin/bash
#SBATCH --ntasks=6
#SBATCH --ntasks-per-socket=6
#SBATCH -t 24:00:00
#
#

export OMP_NUM_THREADS=1

./batch.run.sh 

