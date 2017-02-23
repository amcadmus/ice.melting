#!/bin/bash

nnode=`cat string.out | wc -l`
nnode=$(($nnode-1))

for ii in `seq 0 $nnode`; do lidx=$(($ii+1)); sed -n "${lidx}p" energy.out > `printf %03d $ii`.out; done
