#!/bin/bash

source env.sh
source parameters.sh

cd $system_warming_dir

if echo $system_init_mode | grep traj &> /dev/null; then
    echo "# use conf from traj $system_init_trr at time $system_init_trr_time"
    $exec_grompp -c conf.gro -t $system_init_trr -time $system_init_trr_time &> grompp.log
else 
    if echo $system_init_mode | grep conf &> /dev/null; then
	$exec_grompp -c conf.gro
    fi
fi

if test ! -f topol.tpr; then
    echo "# unsuccessful gromacs preprocess, exit"
    exit
fi

$exec_mdrun -v -plumed plumed.dat

