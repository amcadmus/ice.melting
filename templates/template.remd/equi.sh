#!/bin/bash

source parameters.sh
source env.sh
source $GMX_DIR/bin/GMXRC

dir_base=`pwd`
gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_equi_time) / $gmx_dt" | bc`

nproc=`echo $temperatures | wc -w`

for temperature in $temperatures ;
do
    dir_current=dir.$temperature
    if test ! -d $dir_current; then
	echo "does not exist dir $dir_current, exit"
	exit
    fi
    if test -f $dir_current/equied; then
	echo "# dir $dir_current already equied, continue"
	continue
    fi
    cd $dir_current
    echo "# run in dir $dir_current"
    $grompp_cmd -n index.ndx
    rm -f state.cpt
    $mdrun_equi_cmd -nsteps $gmx_nsteps
    touch equied
    cd $dir_base    
done
    
