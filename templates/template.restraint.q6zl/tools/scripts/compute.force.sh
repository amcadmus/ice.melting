#!/bin/bash

if test $# -lt 1; then
    echo "usage: "
    echo "$0 target_dir [traj_start]"
fi
target_dir=$1
traj_start=10000
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
make -C $compute_force_dir -j4 &>/dev/null

cd $target_dir
if test ! -f parameters.sh; then
    echo "does not exist file parameters.sh, exit"
    exit
fi
source parameters.sh

echo "# compute q6, begin at $traj_start"
command="$compute_force_dir/compute.reweight --input COLVAR.res --column-cv 3 --begin $traj_start --ref-value $plumed_q6_at --evl-value $plumed_q6_at --kappa $plumed_q6_kappa --temperature $gmx_temperature"
echo "# with command: $command"
$command

echo "# compute zl, begin at $traj_start"
command="$compute_force_dir/compute.reweight --input COLVAR.res --column-cv 5 --begin $traj_start --ref-value $plumed_zl_at --evl-value $plumed_zl_at --kappa $plumed_zl_kappa --temperature $gmx_temperature"
echo "# with command: $command"
$command
