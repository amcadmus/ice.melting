#!/bin/bash

# update the ctrl.input and top.input by using parameters.sh

source parameters.sh
source env.sh

if test $# -lt 1 ; then
    echo usage:
    echo $0 dir
    exit 1   
fi
out_dir=$1

tools_dir=`cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd`
cwd=`pwd`

system_start_conf=confs/$system_start_conf
system_ref_conf=$system_start_conf
if test ! -f $system_start_conf; then
    echo no file $system_start_conf
    exit 1
fi

traj_n_freq=`echo "($md_traj_freq + 0.5 * $md_dt) / $md_dt " | bc`
ener_n_freq=`echo "($md_ener_freq + 0.5 * $md_dt) / $md_dt " | bc`
natom=`wc -l $system_ref_conf | awk '{print $1}' `
natom=`echo "($natom - 3)" | bc`
nmol=`echo "$natom / 4" | bc`

rm -f $out_dir/parameters.sh
rm -f $out_dir/top.input
cp parameters.sh $out_dir
top_file=tops/top.input.$system_top_style
cp $top_file $out_dir/top.input

cd $out_dir

# prepare ctrl.input (md parameters)
sed -e "s/\(conf_copies.*=\).*/\1 1,1,1/g" ctrl.input |\
sed -e "s/\(rand_seed.*=\).*/\1 $md_rand_seed/g" |\
sed -e "s/\(end_time.*=\).*/\1 $md_time/g" |\
sed -e "s/\(time_step.*=\).*/\1 $md_dt/g" |\
sed -e "s/\(temperature.*=\).*/\1 $md_T/g" |\
sed -e "s/\(tau_t.*=\).*/\1 $md_tau_t/g" |\
sed -e "s/\(tau_p.*=\).*/\1 $md_tau_p/g" |\
sed -e "s/\(javis_dump_freq.*=\).*/\1 $traj_n_freq/g" |\
sed -e "s/\(energy_freq.*=\).*/\1 $ener_n_freq/g" > tmp.out
mv -f tmp.out ctrl.input

# prepare top.input (restraint parameters)
numb_res=`grep HarmonicRestraint top.input | wc -l`
for ii in `seq 1 $numb_res`
do
    key=`grep HarmonicRestraint top.input | head -n $ii | tail -n 1 | awk '{print $1}'`
    center_value=`echo $res_centers | cut -d ',' -f $ii`
    sed -i "/$key/s/CENTER_VALUE/$center_value/g" top.input
done
sed -i "s/K_VALUE/$res_k/g" top.input
sed -i "/numbers/s/=.*/= $nmol/g" top.input
