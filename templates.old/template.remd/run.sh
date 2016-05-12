#!/bin/bash

source parameters.sh
source env.sh
source $GMX_DIR/bin/GMXRC

dir_base=`pwd`
gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_time) / $gmx_dt" | bc`
gmx_nreplex=`echo "($gmx_dt*0.5+$gmx_replex) / $gmx_dt" | bc`

nproc=`echo $temperatures | wc -w`

for temperature in $temperatures ;
do
    dir_current=dir.$temperature
    if test ! -d $dir_current; then
	echo "does not exist dir $dir_current, exit"
	exit
    fi
    echo "# grompp in dir $dir_current"
    cd $dir_current
    if test ! -f equied; then
	echo "dir $dir_current not equilibriated, exit"
	exit
    fi
    if test ! -f state.cpt; then
	echo "no state.cpt, may not be equilibriated, exit"
	exit
    fi
    $grompp_cmd -n index.ndx -t state.cpt &> grompp.log
    cd $dir_base
done

command="mpirun -n $nproc $mdrun_cmd -multidir dir* -replex $gmx_nreplex -v -nsteps $gmx_nsteps"
echo "#run with command: $command"
$command

