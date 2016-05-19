#!/bin/bash

# compilers
source $HOME/compiler_vars.sh

# batch system
batch_system=PBS
batch_sub=$batch_system.sub

# mpi
mpirun_command=mpirun
num_proc=12

# moasp
moasp_bin_dir=$HOME/tmp/moasp/build/src/app/
moasp_ipp=msp_ipp
moasp_avg=msp_avg
moasp_ctj=msp_cvt_traj
moasp_run=moasp

# gromacs
gmx_bin_dir=$HOME/local/gromacs/4.6.5/bin
gmx_genconf=genconf
gmx_grompp=grompp
gmx_trjconf=trjconv

