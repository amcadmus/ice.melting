#!/bin/bash

source parameters.sh
source env.sh
source $GMX_DIR/bin/GMXRC

dir_base=`pwd`
gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_time) / $gmx_dt" | bc`
gmx_nreplex=`echo "($gmx_dt*0.5+$gmx_replex) / $gmx_dt" | bc`

nproc=`echo $temperatures | wc -w`

command="mpirun -n $nproc $mdrun_cmd -multidir dir* -replex $gmx_nreplex -v -nsteps $gmx_nsteps"
echo "#run with command: $command"
$command

