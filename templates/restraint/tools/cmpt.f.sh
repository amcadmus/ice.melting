#!/bin/bash

source parameters.sh
source env.sh

tools_dir=`cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd`
cwd=`pwd`

numb_res=`grep HarmonicRestraint $top_file | wc -l`
numb_centers=`echo $res_centers | awk -F ',' '{print NF}'`
if test $numb_centers -ne $numb_res; then
    echo "find $numb_res in top while find $numb_centers restraint center, inconsistency, do nothing."
    exit
fi
    
for ii in `seq 1 $numb_res`
do
    key=`grep HarmonicRestraint top.input | head -n $ii | tail -n 1 | awk '{print $1}'`
    col=`grep Print top.input | sed 's/.*input =\(.*\)file.*/\1/g' | sed 's/\"//g' | sed 's/,/\n/g' | grep -n $key | awk -F':' '{print $1}'`
    echo $key $col
done