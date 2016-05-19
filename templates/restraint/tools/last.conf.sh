#!/bin/bash

source env.sh
source $gmx_bin_dir/GMXRC.bash

target_dir="."
if [ $# -ge 1 ]; then
    target_dir=$1
fi
if [ ! -d $target_dir ]; then
    echo "no target $target_dir, exit"
    exit
fi
viz_dir=$target_dir/viz
if [ ! -d $viz_dir ]; then
    echo "no viz $viz_dir, exit"
    exit
fi

frame_id=`tail $viz_dir/REAL_SPACE/dumps.javis -n 1 | cut -d '/' -f 1 | cut -d '.' -f 2`

cd $target_dir
$moasp_bin_dir/$moasp_ctj -n $frame_id -m trr
if [ $? -ne 0 ]; then echo "unsuccessful msp_cvt_traj, exit"; exit; fi

cd gmx.tool
$gmx_grompp  &> /dev/null
if [ $? -ne 0 ]; then echo "unsuccessful grompp, exit"; exit; fi

echo 0 | $gmx_trjconf -f ../out.trr -o ../out.gro -pbc mol &> /dev/null
if [ $? -ne 0 ]; then echo "unsuccessful trajconv, exit"; exit; fi

cd ..
