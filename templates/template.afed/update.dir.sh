#!/bin/bash

source parameters.sh

gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_time) / $gmx_dt" | bc`
gmx_nxtc=`echo "($gmx_dt*0.5+$gmx_xtc_feq) / $gmx_dt" | bc`
gmx_ntrr=`echo "($gmx_dt*0.5+$gmx_trr_feq) / $gmx_dt" | bc`
gmx_nenergy=`echo "($gmx_dt*0.5+$gmx_energy_feq) / $gmx_dt" | bc`

dir_base=`pwd`

dir_current=$system_running_dir

if test ! -d $dir_current; then
    echo "does not exist dir $dir_current, exit"
    exit
fi

echo "# update parameters in $dir_current"
    
cd $dir_current
tau_p_line="$gmx_tau_box $gmx_tau_box $gmx_tau_box $gmx_tau_box $gmx_tau_box $gmx_tau_box"
ref_p_line="1 1 1 0 0 0"
compres_line="4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5"
sed -e "s/dt.*=.*/dt = $gmx_dt/g" grompp.mdp | \
    sed -e "s/nstlog.*=.*/nstlog = 0/g" | \
    sed -e "s/nstenergy.*=.*/nstenergy = $gmx_nenergy/g" | \
    sed -e "s/nstcomm.*=.*/nstcomm = $gmx_nenergy/g" | \
    sed -e "s/nstxtcout.*=.*/nstxtcout = $gmx_nxtc/g" | \
    sed -e "s/nstxout.*=.*/nstxout = $gmx_ntrr/g" | \
    sed -e "s/nstvout.*=.*/nstvout = $gmx_ntrr/g" | \
    sed -e "s/nstfout.*=.*/nstfout = 0/g" | \
    sed -e "s/gen_vel.*=.*/gen_vel = no/g" | \
    sed -e "s/ref_t.*=.*/ref_t = $gmx_temperature/g" |\
    sed -e "s/tau_t.*=.*/tau_t = $gmx_tau_t/g" |\
    sed -e "s/Pcoupltype.*=.*/Pcoupltype = anisotropic/g" |\
    sed -e "s/ref_p.*=.*/ref_p = $ref_p_line/g" |\
    sed -e "s/tau_p.*=.*/tau_p = $tau_p_line/g" |\
    sed -e "s/compressibility.*=.*/compressibility = $compres_line/g" |\
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


