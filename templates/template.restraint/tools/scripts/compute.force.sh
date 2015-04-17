#!/bin/bash

if test $# -lt 1; then
    echo "usage: "
    echo "$0 target_dir [traj_start]"
fi
target_dir=$1
traj_start=500
if test $# -ge 2; then
    traj_start=$2
fi

if test ! -d $target_dir; then
    echo "does not exist dir $target_dir, exit"
    exit
fi
echo "# compute the force in dir $target_dir"

script_dir=$(cd ${0%/*} && echo $PWD)
base_dir=$script_dir/../
compute_force_dir=$base_dir/compute.force
echo $script_dir
echo $compute_force_dir

cd $target_dir
if test ! -f parameters.sh; then
    echo "does not exist file parameters.sh, exit"
    exit
fi
source parameters.sh
echo "# compute q4, begin at $traj_start"
$compute_force_dir --input COVAR.res --column-cv 2 --begin $traj_start --ref-value $plumed_q4_at --evl-value $plumed_q4_at --kappa $plumed_q4_kappa --temperature $gmx_temperature
echo "# compute q6, begin at $traj_start"
$compute_force_dir --input COVAR.res --column-cv 3 --begin $traj_start --ref-value $plumed_q6_at --evl-value $plumed_q6_at --kappa $plumed_q6_kappa --temperature $gmx_temperature
