#!/bin/bash

source parameters.sh

gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_time) / $gmx_dt" | bc`
gmx_nxtc=`echo "($gmx_dt*0.5+$gmx_xtc_feq) / $gmx_dt" | bc`
gmx_ntrr=`echo "($gmx_dt*0.5+$gmx_trr_feq) / $gmx_dt" | bc`
gmx_nenergy=`echo "($gmx_dt*0.5+$gmx_energy_feq) / $gmx_dt" | bc`

dir_base=`pwd`

for temperature in $temperatures ;
do 
    dir_current=dir.$temperature
    if test -d $dir_current; then
	echo "exist dir $dir_current, continue"
	continue
    fi
    if test ! -d md.seed; then
	echo "the dir md.seed does not exist, exit"
	exit
    fi

    echo "# generate dir $dir_current"
    
    # echo "## generate mdp"
    cp -a md.seed $dir_current    

    cd $dir_base
done

./update.dir.sh

