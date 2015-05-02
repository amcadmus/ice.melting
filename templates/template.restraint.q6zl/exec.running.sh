#!/bin/bash

source env.sh
source parameters.sh

rm -f $system_running_dir/conf.gro
cp $system_warming_dir/confout.gro  $system_running_dir/conf.gro

cd $system_running_dir/

$exec_grompp -c conf.gro

if test ! -f topol.tpr; then
    echo "# unsuccessful gromacs preprocess, exit"
    exit
fi

$exec_mdrun -v -plumed plumed.dat

