#!/bin/bash

moasp_dir=~/study/moasp/build/src/app/

if test $# -lt 1 ; then
    echo "# print the histogram of mol value "
    echo "# usage"
    echo $0 dir [dir [dir ...] ]
    exit
fi

targets=$*

for dir in $targets
do
    if test ! -d $dir; then
	echo no dir $dir
	exit
    fi

    cwd=`pwd`
    cd $dir/
    cat mol_* | grep -v \# > all.mol.out;
    $moasp_dir/msp_histo -f all.mol.out -u .7 -l .0 -b 200 -o hist.mol.out -c 3; 
    cd $cwd
done

