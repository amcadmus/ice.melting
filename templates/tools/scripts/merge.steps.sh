#!/bin/bash

if test $# -ne 3; then
    echo usage
    echo $0 input0 input1 output
    exit
fi

in0=$1
in1=$2
output=$3

if test ! -d $in0; then
    echo nodir $in0; 
    exit
fi
if test ! -d $in1; then
    echo nodir $in1; 
    exit
fi

test ! -d $output && mkdir -p $output

cwd=`pwd`

echo input0: $in0
echo input1: $in1
echo output: $output
cd $in0
targets=`ls | grep ^step | grep [0-9]$`
cd $cwd
echo will merge files: $targets

for ii in $targets; 
do
    join $in0/$ii $in1/$ii > $output/$ii
done


