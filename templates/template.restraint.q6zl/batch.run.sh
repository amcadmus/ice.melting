#!/bin/bash

if test $# -lt 3; then
    echo "usage:"
    echo "$0 q6 zl start_time"
fi

q6=$1
zl=$2
start_time=$3
q6_ptr=`printf %.3f $q6`
zl_ptr=`printf %.3f $zl`
warming_dir=warming.$q6_ptr.$zl_ptr
running_dir=running.$q6_ptr.$zl_ptr

sed -e "s/plumed_q6_at=.*/plumed_q6_at=$q6/g" parameters.sh |\
sed -e "s/plumed_zl_at=.*/plumed_zl_at=$zl/g" |\
sed -e "s/system_warming_dir=.*/system_warming_dir=$warming_dir/g" |\
sed -e "s/system_running_dir=.*/system_running_dir=$running_dir/g" |\
sed -e "s/system_init_trr_time=.*/system_init_trr_time=$start_time/g" > tmp.tmp
mv -f tmp.tmp parameters.sh

./gen.dir.sh
./exec.warming.sh
./exec.running.sh


