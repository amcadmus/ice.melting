#!/bin/bash

source env.sh
source parameters.sh

source $GMX_DIR/bin/GMXRC

cd $system_warming_dir

if echo $system_init_mode | grep traj &> /dev/null; then
    echo "# use conf from traj $system_init_xtc at time $system_init_xtc_time"
    trjconv -b $system_init_xtc_time -e $system_init_xtc_time -o conf.gro
    grompp_mpi -c conf.gro
else 
    if echo $system_init_mode | grep conf &> /dev/null; then
	grompp_mpi -c conf.gro
    fi
fi

if test ! -f topol.tpr; then
    echo "# unsuccessful gromacs preprocess, exit"
    exit
fi

mpirun -n 2 mdrun_mpi -v -plumed plumed.dat

