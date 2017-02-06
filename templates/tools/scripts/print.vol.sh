#!/bin/bash

moasp_dir=~/study/moasp/build/src/app/

if test $# -lt 1 ; then
    echo "# print volume of the string node "
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
    result=`$moasp_dir/msp_avg -f $dir/energy.out -c 9 -t .8 | grep -v \#`
    echo $dir $result
done

