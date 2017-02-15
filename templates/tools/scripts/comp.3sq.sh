#!/bin/bash

if test $# -lt 1; then
    echo usage
    echo $0 [ string [ string [ string ... ] ] ]
    exit
fi

targets=$*
cwd=`pwd`
sq_dir=/home/wanghan/study/ice.melting/templates/tools/sq

for ii in $targets
do
    if test ! -d $ii; then continue; fi
    cd $ii;
    for jj in node.0*; 
    do
	s100=`$sq_dir/sq.single -f $jj/out.xtc --qx 0 --qy 16.158 --qz 0 | grep -v \# `
	s002=`$sq_dir/sq.single -f $jj/out.xtc --qx 0 --qy 0 --qz 17.217 | grep -v \#`
	s101=`$sq_dir/sq.single -f $jj/out.xtc --qx 0 --qy 16.158 --qz 8.609 | grep -v \#`
	echo $s100 $s002 $s101
    done
    cd $cwd
done
