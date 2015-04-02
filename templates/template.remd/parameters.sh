#!/bin/bash

gmx_equi_time=10		# in ps
gmx_time=1000000		# in ps
gmx_replex=1			# in ps

gmx_dt=0.002			# in ps
gmx_xtc_feq=10.0		# in ps
gmx_trr_feq=0			# in ps

grompp_cmd=grompp
mdrun_equi_cmd="mdrun -nt 4 -v"
mdrun_cmd="mdrun_mpi"
#gmx_dir=$HOME/local/gromacs/4.6.5

temperatures="150.00 155.11 160.34 165.67 171.13 176.71 182.40 188.22 194.17 200.24 206.45 212.80 219.29 225.93 232.72 239.66 246.75 252.10 253.88 261.28 268.85 276.60 284.53 292.64 300.92 309.40 318.06 326.93 336.00 345.28 354.76"

