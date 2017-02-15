#!/bin/bash

target_dir="."
if [ $# -ge 1 ]; then
    target_dir=$1
fi
if [ ! -d $target_dir ]; then
    echo "no target $target_dir, exit"
    exit 1
fi

tools_dir=`cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd`
cwd=`pwd`

cd $target_dir
source env.sh
source parameters.sh

numb_res=`grep "HarmonicRestraint\|VolumeConstraint" top.input | wc -l`
numb_centers=`echo $res_centers | awk -F ',' '{print NF}'`
if test $numb_centers -ne $numb_res; then
    echo "find $numb_res in top while find $numb_centers restraint center, inconsistency, do nothing."
    exit 1
fi

line=""
for ii in `seq 1 $numb_res`
do
    param_line=`grep "HarmonicRestraint\|VolumeConstraint" top.input  | head -n $ii | tail -n 1`
    key=`echo $param_line | sed 's/.*input =\(.*\)k.*/\1/g' | awk -F'"' '{print $2}'`
    center=`echo $param_line | sed 's/.*center =\(.*\)bias.*/\1/g'`
    kk=`echo $param_line | sed 's/.*k =\(.*\)center.*/\1/g'`
    col=`sed 's;//.*;;g' top.input | grep Print | sed 's/.*input =\(.*\)file.*/\1/g' | sed 's/\"//g' | sed 's/,/\n/g' | grep -n $key | awk -F ':' '{print $1}'`
    col=`echo "$col+1" | bc`   
    echo "# $key $col $center $kk"
    awk "{print (\$$col - $center) * $kk}" mace.out > tmp.out
    $moasp_bin_dir/$moasp_avg -f tmp.out -t .8 > avg.out
    force=`grep -v \# avg.out | awk '{print $1}'`
    force_err=`grep -v \# avg.out | awk '{print $2}'`
    if [ $? -ne 0 ]; then echo "unsuccessful msp_avg, exit"; exit 1; fi
    line="$line $force"
    line_err="$line_err $force_err"
    line_center="$line_center $center"
done
rm -f tmp.out

cd $cwd
echo $line
echo "$line $line_err" > node.force.out
echo "$line_center $line $line_err" > node.all.out
