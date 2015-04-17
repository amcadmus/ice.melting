#!/bin/bash

if test $# -lt 3; then
    echo "usage:"
    echo "$0 q4 q6 start_time"
fi

q4=$1
q6=$2
start_time=$3
q4_ptr=`printf %.3f $q4`
q6_ptr=`printf %.3f $q6`
warming_dir=warming.$q4_ptr.$q6_ptr
running_dir=running.$q4_ptr.$q6_ptr

sed -e 's/plumed_q4_at=.*/plumed_q4_at=$q4/g' parameters.sh |\
sed -e 's/plumed_q6_at=.*/plumed_q6_at=$q6/g' |\
sed -e 's/system_warming_dir=.*/system_warming_dir=$warming_dir/g' |\
sed -e 's/system_running_dir=.*/system_running_dir=$running_dir/g' |\
sed -e 's/system_init_xtc_time=.*/system_init_xtc_time=$start_time/g' > tmp.tmp
mv -f tmp.tmp parameters.sh

./gen.dir.sh
./exec.warming.sh
./exec.running.sh


