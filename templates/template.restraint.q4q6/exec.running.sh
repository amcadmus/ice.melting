#!/bin/bash

source env.sh
source parameters.sh

source $GMX_DIR/bin/GMXRC

rm -f $system_running_dir/conf.gro
cp $system_warming_dir/confout.gro  $system_running_dir/conf.gro

cd $system_running_dir/

grompp_mpi -c conf.gro

if test ! -f topol.tpr; then
    echo "# unsuccessful gromacs preprocess, exit"
    exit
fi

mpirun -n 2 mdrun_mpi -v -plumed plumed.dat

