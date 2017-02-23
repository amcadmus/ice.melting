#!/bin/bash

moasp_dir=~/study/moasp/build/src/app/
cv_dir=/home/wanghan/study/ice.melting/templates/tools/cv

if test $# -lt 1 ; then
    echo "# print volume of the string node "
    echo "# usage"
    echo $0 dir [dir [dir ...] ]
    exit
fi

targets=$*

cwd=`pwd`
for dir in $targets
do
    if test ! -d $dir; then
	echo no dir $dir
	exit
    fi
    if test ! -f $dir/out.xtc ; then
	echo no file $dir/out.xtc, skip
	continue
    fi
    cd $dir
    result=`$cv_dir/locvol -f out.xtc | grep -v \#`
    cd $cwd
    echo $dir $result
done

