#!/bin/bash

target_dir="."
if [ $# -ge 1 ]; then
    target_dir=$1
fi
if [ ! -d $target_dir ]; then
    echo "no target $target_dir, exit"
    exit 1
fi
viz_dir=$target_dir/viz
if [ ! -d $viz_dir ]; then
    echo "no viz $viz_dir, exit"
    exit 1
fi

cd $target_dir
source env.sh
source $gmx_bin_dir/GMXRC.bash

frame_id=`tail $viz_dir/REAL_SPACE/dumps.javis -n 1 | cut -d '/' -f 1 | cut -d '.' -f 2`

$moasp_bin_dir/$moasp_ctj -n $frame_id -m trr
if [ $? -ne 0 ]; then echo "unsuccessful msp_cvt_traj, exit"; exit 1; fi

cd gmx.tool
# $gmx_grompp  &> /dev/null
# if [ $? -ne 0 ]; then echo "unsuccessful grompp, exit"; exit 1; fi

echo 0 | $gmx_trjconf -f ../out.trr -o ../out.gro -pbc mol &> /dev/null
if [ $? -ne 0 ]; then echo "unsuccessful trajconv, exit"; exit 1; fi

cd ..
