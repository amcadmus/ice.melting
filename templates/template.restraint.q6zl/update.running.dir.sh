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

echo "# update grompp.mdp"
if echo $gmx_pcoupltype | grep anisotropic &> /dev/null; then
    tau_p_line="$gmx_tau_box $gmx_tau_box $gmx_tau_box $gmx_tau_box $gmx_tau_box $gmx_tau_box"
    ref_p_line="1 1 1 0 0 0"
    compres_line="4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5 4.5e-5"
else 
    tau_p_line="$gmx_tau_box"
    ref_p_line="1"
    compres_line="4.5e-5"
fi
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
    sed -e "s/Pcoupltype .*=.*/Pcoupltype = $gmx_pcoupltype/g" |\
    sed -e "s/ref_p.*=.*/ref_p = $ref_p_line/g" |\
    sed -e "s/tau_p.*=.*/tau_p = $tau_p_line/g" |\
    sed -e "s/compressibility.*=.*/compressibility = $compres_line/g" |\
    sed -e "s/nsteps.*=.*/nsteps = $gmx_nsteps/g" > tmp.tmp
mv -f tmp.tmp grompp.mdp    

echo "# update plumed.dat"
rm -f plumed.dat
cp $dir_base/md.seed/plumed.dat .
nline=`wc -l conf.gro | awk '{print $1}'`
nline=$(($nline-3))
nmol=`echo "$nline / 4" | bc`
ratio=1
for ii in `echo $system_box`;
do
    ratio=`echo "$ratio * $ii" | bc`
done
effective_q6_kappa=`echo "$ratio * $plumed_q6_kappa" | bc -l`
effective_zl_kappa=`echo "$ratio * $plumed_zl_kappa" | bc -l`
nbox_z=`echo $system_box | awk '{print $3}' `
effective_n_zl=`echo "$plumed_zl_ref_numb * $nbox_z" | bc `
sed -e "s/INPUT_NATOM/$nline/g" plumed.dat |\
    sed -e "s/INPUT_NLIST_RC/$plumed_nlist_rc/g" |\
    sed -e "s/INPUT_NLIST_FEQ/$plumed_nlist_feq/g" |\
    sed -e "s/INPUT_Q6/$plumed_q6_at/g" |\
    sed -e "s/INPUT_ZL/$plumed_zl_at/g" |\
    sed -e "s/INPUT_KAPPA_Q6/$effective_q6_kappa/g" |\
    sed -e "s/INPUT_KAPPA_ZL/$effective_zl_kappa/g" |\
    sed -e "s/INPUT_NUMB_ZL/$effective_n_zl/g" |\
    sed -e "s/INPUT_COMP_SHIFT/$plumed_zl_comp_shift/g" |\
    sed -e "s/INPUT_DT/$gmx_warm_dt/g" > tmp.tmp
mv -f tmp.tmp plumed.dat

# echo "## generate index"
echo "# update index.ndx"
rm -f index.ndx
if test ! -d index.ndx; then
    echo "a OW" > tmp.in
    echo "q" >> tmp.in
    cat tmp.in | make_ndx -f  conf.gro &> /dev/null
    rm -f tmp.in
fi

echo "# backup parameter"
rm -f parameters.sh
cp $dir_base/parameters.sh .

cd $dir_base


