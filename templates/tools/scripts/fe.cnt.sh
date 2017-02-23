#!/bin/bash

surface_ener=27
defect_file=steinhardt.avg.def.out
moasp_bin=$HOME/local/bin/

if test $# -eq 0; then
    echo usage
    echo $0 [string [string [string ...] ] ]
    exit
fi

targets=$*

cwd=`pwd`
for ii in $targets
do
    cd $ii
    nnode=`cat string.out | wc -l`
    nnode=$(($nnode-1))
    
    rm -f fe.cnt.out
    for jj in `seq 0 $nnode`
    do
	node_dir=node.`printf %06d $jj`
	if test ! -f $node_dir/$defect_file; then
	    continue
	fi
	numb_def=`$moasp_bin/msp_avg -f $node_dir/$defect_file -c 3 | grep -v \# | awk '{print $1}'`
	numb_mol=`grep numbers $node_dir/top.input | grep = | cut -d '=' -f 2`
	vol=`$moasp_bin/msp_avg -f $node_dir/energy.out -c 9 | grep -v \# | awk '{print $1}'`
	def_vol=`echo "$numb_def $numb_mol $vol" | awk '{print $1/$2*$3}'`
	def_radius=`echo "$def_vol" | awk '{print ($1 / (4./3.*3.1416))^(1./3.)}'`
	def_surf=`echo $def_radius | awk '{print 4.*3.1416*$1*$1}'`
#	echo $def_surf
	def_ener=`echo "$def_surf $surface_ener" | awk '{print $1*$2/1.661}'`
	echo $jj $def_ener >> fe.cnt.out
    done
    cd $cwd
done


	      
	      
	      

