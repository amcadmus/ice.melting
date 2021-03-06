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
    if test ! -d $dir_current; then
	echo "does not exist dir $dir_current, exit"
	exit
    fi

    echo "# update dir $dir_current"
    
    # echo "## generate mdp"
    cd $dir_current
    sed -e "s/ref_t.*=.*/ref_t = $temperature/g" grompp.mdp | \
    sed -e "s/nstlog.*=.*/nstlog = 0/g" | \
    sed -e "s/nstenergy.*=.*/nstenergy = $gmx_nenergy/g" | \
    sed -e "s/nstcomm.*=.*/nstcomm = $gmx_nenergy/g" | \
    sed -e "s/nstxtcout.*=.*/nstxtcout = $gmx_nxtc/g" | \
    sed -e "s/nstxout.*=.*/nstxout = $gmx_ntrr/g" | \
    sed -e "s/nstvout.*=.*/nstvout = $gmx_ntrr/g" | \
    sed -e "s/nstfout.*=.*/nstfout = 0/g" | \
    sed -e "s/gen_vel.*=.*/gen_vel = no/g" | \
    sed -e "s/nsteps.*=.*/nsteps = $gmx_nsteps/g" > tmp.tmp
    mv -f tmp.tmp grompp.mdp    

#    sed -e "s/xtc-grps.*=.*/xtc-grps = OW/g" | \

    # echo "## generate index"
    if test ! -d index.ndx; then
	echo "a OW" > tmp.in
	echo "q" >> tmp.in
	cat tmp.in | make_ndx -f  conf.gro &> /dev/null
	rm -f tmp.in
    fi

    cd $dir_base
done


