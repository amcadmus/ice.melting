#!/bin/bash

HEAD="==run_info=="
# check if finished
if [ -f tag_finished ]; then
    echo "$HEAD finished job, do nothing"
    echo "$HEAD if want to start anyway, remove tag: tag_finished"
    exit 1
fi

source env.sh
export OMP_NUM_THREADS=1

# ipp
command="$moasp_bin_dir/$moasp_ipp -n $numb_proc"
echo "$HEAD ipp with command $command"
$command
if [ $? -ne 0 ]; then
    echo "failed to tun $moasp_ipp, exit"
    exit 1
fi

# md run
command="$mpirun_command -n $numb_proc $moasp_bin_dir/$moasp_run -v &> /dev/null"
echo "$HEAD run with command $command"
$command
if [ $? -eq 0 ]; then
    echo "$HEAD finishes successfully"
    touch tag_finished
else
    echo "$HEAD finishes unsuccessfully"
fi
