#!/bin/bash

source env.sh
source parameters.sh

source $GMX_DIR/bin/GMXRC

cd $system_warming_dir

if echo $system_init_mode | grep traj &> /dev/null; then
    echo "# use conf from traj $system_init_trr at time $system_init_trr_time"
    grompp_mpi -c conf.gro -t $system_init_trr -time $system_init_trr_time &> grompp.log
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

