#!/bin/bash

HEAD="==run_info=="
# check if finished
if [ -f tag_finished ]; then
    echo "$HEAD finished job, do nothing"
    echo "$HEAD if want to start anyway, remove tag: tag_finished"
    exit
fi

source env.sh
export OMP_NUM_THREADS=1

# ipp
command="$moasp_bin_dir/$moasp_ipp -n $num_proc"
echo "$HEAD ipp with command $command"
$command
if [ $? -ne 0 ]; then
    echo "failed to tun $moasp_ipp, exit"
    exit
fi

# md run
command="$mpirun_command -n $num_proc $moasp_bin_dir/$moasp_run -v &> /dev/null"
echo "$HEAD run with command $command"
$command
if [ $? -eq 0 ]; then
    echo "$HEAD finishes successfully"
    touch tag_finished
else
    echo "$HEAD finishes unsuccessfully"
fi
