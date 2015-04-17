#!/bin/bash

source env.sh
source parameters.sh

source $GMX_DIR/bin/GMXRC

cd $system_warming_dir

if echo $system_init_mode | grep traj &> /dev/null; then
    echo "# use conf from traj $system_init_xtc at time $system_init_xtc_time"
    grompp_mpi -c conf.gro &> grompp.0.log
    echo 0 | trjconv -f $system_init_xtc -b $system_init_xtc_time -e $system_init_xtc_time -o conf.gro &> trjconv.log
    grompp_mpi -c conf.gro &> grompp.log
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

