#!/bin/bash

source env.sh

export OMP_NUM_THREADS=1

$moasp_bin_dir/msp_ipp -n $num_proc

command="$mpirun_command -n $num_proc $moasp_bin_dir/moasp -v &> /dev/null"

echo "run with command $command"

$command
