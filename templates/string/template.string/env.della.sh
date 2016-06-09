#!/bin/bash

# compilers
source $HOME/module.sh

# batch system
batch_system=Slurm
batch_sub=$batch_system.sub
batch_queue=
job_hour=1
job_min=0
cput_hour=

# mpi
mpirun_command=srun
numb_proc=15
numb_proc_per_node=5

# moasp
moasp_bin_dir=$HOME/local/moasp/bin/
moasp_ipp=msp_ipp
moasp_avg=msp_avg
moasp_ctj=msp_cvt_traj
moasp_run=moasp

# gromacs
gmx_bin_dir=$HOME/local/gromacs/4.6.7/bin
gmx_genconf=genconf
gmx_grompp=grompp
gmx_trjconf=trjconv

