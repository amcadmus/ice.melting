#!/bin/bash

moasp_dir=~/study/moasp/build/src/app/

if test $# -lt 2 ; then
    echo "# smooths mole traj files matches input/mol_* and dumps to output/mol_*"
    echo "# usage"
    echo $0 input output [support] [column]
    echo "# with default support 5"
    echo "# with default column  3"
    exit
fi

support=5
col=3
if test $# -ge 2; then
    input=$1
    output=$2
fi
if test $# -ge 3; then
    support=$3
fi
if test $# -ge 4; then
    col=$4
fi

echo "# in: $input"
echo "# out: $output"
echo "# support: $support"

if test ! -d $output; then
    mkdir -p $output
fi

target=$input/mol_*
for ii in $target;
do
    ofile=`echo $ii | sed "s;$input;$output;g"`
    printf "# convert $ii to $ofile \r"
    $moasp_dir/msp_smth -f $ii -m $col -r $support > $ofile
done
echo ""

