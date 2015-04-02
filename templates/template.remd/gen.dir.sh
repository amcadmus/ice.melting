#!/bin/bash

source parameters.sh

gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_time) / $gmx_dt" | bc`
gmx_nxtc=`echo "($gmx_dt*0.5+$gmx_xtc_feq) / $gmx_dt" | bc`
gmx_ntrr=`echo "($gmx_dt*0.5+$gmx_trr_feq) / $gmx_dt" | bc`

dir_base=`pwd`

for temperature in $temperatures ;
do 
    dir_current=dir.$temperature
    if test -d $dir_current; then
	echo "exist dir $dir_current, exit"
	exit
    fi
    if test ! -d md.seed; then
	echo "the dir md.seed does not exist, exit"
	exit
    fi

    echo "# generate dir $dir_current"
    
    # echo "## generate mdp"
    cp -a md.seed $dir_current    
    cd $dir_current
    sed -e "s/ref_t.*=.*/ref_t = $temperature/g" grompp.mdp | \
    sed -e "s/xtc-grps.*=.*/xtc-grps = OW/g" | \
    sed -e "s/nstxtcout.*=.*/nstxtcout = $gmx_nxtc/g" | \
    sed -e "s/nstxout.*=.*/nstxout = $gmx_ntrr/g" | \
    sed -e "s/nstvout.*=.*/nstvout = $gmx_ntrr/g" | \
    sed -e "s/nstfout.*=.*/nstfout = 0/g" | \
    sed -e "s/gen_vel.*=.*/gen_vel = no/g" | \
    sed -e "s/nsteps.*=.*/nsteps = $gmx_nsteps/g" > tmp.tmp
    mv -f tmp.tmp grompp.mdp    

    # echo "## generate index"
    if test ! -d index.ndx; then
	echo "a OW" > tmp.in
	echo "q" >> tmp.in
	cat tmp.in | make_ndx -f  conf.gro &> /dev/null
	rm -f tmp.in
    fi

    cd $dir_base
done


