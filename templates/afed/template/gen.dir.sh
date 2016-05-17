#!/bin/bash

source parameters.sh
source $gmx_bin_dir/GMXRC.bash

if test ! -d $system_seed_dir; then
    echo "no dir $system_seed_dir"
    exit
fi
if test ! -f tops/top.input.$system_top_style; then
    echo "no top file tops/top.input.$system_top_style, do nothing."
    top_files=`ls tops | grep ^top`
    echo "available top files are"
    echo "$top_files"
    exit
fi

copy_x=`echo $system_conf_copies | cut -d ',' -f 1`
copy_y=`echo $system_conf_copies | cut -d ',' -f 2`
copy_z=`echo $system_conf_copies | cut -d ',' -f 3`
copy_x=`printf %02d $copy_x`
copy_y=`printf %02d $copy_y`
copy_z=`printf %02d $copy_z`
afed_tau=`printf %04d $afed_tau`
out_dir=${system_start_state}.${copy_x}x${copy_y}x${copy_z}.t${afed_T}K.tau$afed_tau

if test -d $out_dir; then
    echo "existing dir $out_dir, do nothing"
    exit
fi

traj_n_freq=`echo "($md_traj_freq + 0.5 * $md_dt) / $md_dt " | bc`
ener_n_freq=`echo "($md_ener_freq + 0.5 * $md_dt) / $md_dt " | bc`

# copy dir
echo "# using conf file md/$system_start_conf"
echo "# using top file tops/top.input.$system_top_style"
cp -a $system_seed_dir $out_dir
rm -f $out_dir/top.input
cp tops/top.input.$system_top_style $out_dir/top.input
cp parameters.sh $out_dir
cp run.sh env.sh $out_dir
cd $out_dir

# link conf
ln -s $system_start_conf conf.gro

# change ctrl.input (md parameters)
sed -e "s/\(conf_copies.*=\).*/\1 $copy_x,$copy_y,$copy_z/g" ctrl.input |\
sed -e "s/\(rand_seed.*=\).*/\1 $md_rand_seed/g" |\
sed -e "s/\(end_time.*=\).*/\1 $md_time/g" |\
sed -e "s/\(time_step.*=\).*/\1 $md_dt/g" |\
sed -e "s/\(temperature.*=\).*/\1 $md_T/g" |\
sed -e "s/\(tau_t.*=\).*/\1 $md_tau_t/g" |\
sed -e "s/\(tau_p.*=\).*/\1 $md_tau_p/g" |\
sed -e "s/\(javis_dump_freq.*=\).*/\1 $traj_n_freq/g" |\
sed -e "s/\(energy_freq.*=\).*/\1 $ener_n_freq/g" > tmp.out
mv -f tmp.out ctrl.input

# change top.input (afed parameters)
sed -e "/AFED/s/temperature =.*tau/temperature = $afed_T tau/g" top.input |\
sed -e "/AFED/s/tau_t = [0-9]*/tau_t = $afed_tau/g" > tmp.out
mv -f tmp.out top.input

# update gromacs tools
cd tool.gmx
$gmx_genconf -nbox $copy_x $copy_y $copy_z -f ../$system_start_conf -o conf.gro &> /dev/null
gmx_natom=`grep SOL topol.top | awk '{print $2}'`
gmx_natom=`echo "$gmx_natom * $copy_x * $copy_y * $copy_z" | bc`
sed -i "s/SOL.*[0-9]*/SOL $gmx_natom/g" topol.top
$gmx_grompp &> /dev/null




    

