#!/bin/bash

source parameters.sh

dir_base=`pwd`
gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_equi_time) / $gmx_dt" | bc`

for temperature in $temperatures ;
do 
    dir_current=dir.$temperature
    if test ! -d $dir_current; then
	echo "does not exist dir $dir_current, exit"
	exit
    fi
    
    echo "# run in dir $dir_current"
    cd $dir_current
    $grompp_cmd -n index.ndx
    $mdrun_equi_cmd -nsteps $gmx_nsteps
    cd $dir_base
done
    
