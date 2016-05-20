#!/bin/bash

source parameters.sh
source env.sh
source $gmx_bin_dir/GMXRC.bash

tools_dir=`cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd`
cwd=`pwd`
log=gen.dir.log
rm -f $log

if test ! -d $system_seed_dir; then
    echo "no dir $system_seed_dir"
    exit 1
fi
if test ! -f tops/top.input.$system_top_style; then
    echo "no top file tops/top.input.$system_top_style, do nothing."
    top_files=`ls tops | grep ^top`
    echo "available top files are"
    echo "$top_files"
    exit 1
fi
if test $# -ge 1; then
    init_folder=$1
    if test ! -d $init_folder; then
	echo "no init folder $init_folder, do nothing"
	exit 1
    fi
    if test ! -f $init_folder/tag_finished; then
	echo "simulation of $init_folder/ is not finished, do nothing"
	exit 1
    fi
    if test ! -f $init_folder/out.gro; then
	echo "no last conf $init_folder/out.gro, try to generate one"
	cd $init_folder
	$tools_dir/last.conf.sh
	cd $cwd
    fi
    cd $init_folder
    init_folder=`pwd`
    cd $cwd
    system_start_conf=$init_folder/out.gro	
    system_conf_copies=1,1,1			# ncopy reset if start from finished trajs
    revise_top=1
else
    init_folder=$system_seed_dir
    cd $init_folder
    init_folder=`pwd`
    cd $cwd
    system_start_conf=$init_folder/conf.${system_start_state}.gro	# in $system_seed_dir
    revise_top=0
fi

# copy and output freq
copy_x=`echo $system_conf_copies | cut -d ',' -f 1`
copy_y=`echo $system_conf_copies | cut -d ',' -f 2`
copy_z=`echo $system_conf_copies | cut -d ',' -f 3`
copy_x=`printf %02d $copy_x`
copy_y=`printf %02d $copy_y`
copy_z=`printf %02d $copy_z`
natom=`wc -l $system_start_conf | awk '{print $1}' `
natom=`echo "($natom - 3) / 4 * $copy_x * $copy_y * $copy_z" | bc`
natom=`printf %07d $natom`
traj_n_freq=`echo "($md_traj_freq + 0.5 * $md_dt) / $md_dt " | bc`
ener_n_freq=`echo "($md_ener_freq + 0.5 * $md_dt) / $md_dt " | bc`

# parse constraints.
top_file=tops/top.input.$system_top_style
numb_res=`grep HarmonicRestraint $top_file | wc -l`
numb_centers=`echo $res_centers | awk -F ',' '{print NF}'`
if test $numb_centers -ne $numb_res; then
    echo "find $numb_res in top while find $numb_centers restraint center, inconsistency, do nothing."
    cd ..
    rm -fr $out_dir
    exit 1
fi
res_name=""
for ii in `seq 1 $numb_centers`
do
    key=`grep HarmonicRestraint $top_file | head -n $ii | tail -n 1 | awk '{print $1}'`
    center_value=`echo $res_centers | cut -d ',' -f $ii`
    center_value=`printf %.04f $center_value`
    res_name=$res_name.${key}$'-'${center_value}
done
out_dir=${dir_prefix}${system_start_state}.${natom}$res_name

# make output dir
if test -d $out_dir; then
    echo "existing dir $out_dir, do nothing"
    exit 1
fi
# copy dir
echo "# generate dir $out_dir"			>> $log
echo "# using conf file $system_start_conf"	>> $log
echo "# using top file $top_file"		>> $log
cp -a $system_seed_dir $out_dir
rm -f $out_dir/conf.*.gro
rm -f $out_dir/top.input
cp $top_file $out_dir/top.input
cp parameters.sh $out_dir
cp env.sh $out_dir
cd $out_dir

# prepare conf
ln -s $system_start_conf conf.gro

# prepare ctrl.input (md parameters)
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

# prepare top.input (restraint parameters)
for ii in `seq 1 $numb_res`
do
    key=`grep HarmonicRestraint top.input | head -n $ii | tail -n 1 | awk '{print $1}'`
    center_value=`echo $res_centers | cut -d ',' -f $ii`
    sed -i "/$key/s/CENTER_VALUE/$center_value/g" top.input
done
sed -i "s/K_VALUE/$res_k/g" top.input
if [ $revise_top -eq 1 ] ; then
    echo $natom
    sed -i "/numbers/s/=.*/= $natom/g" top.input
fi

# update gromacs tools
cd gmx.tool
$gmx_genconf -nbox $copy_x $copy_y $copy_z -f $system_start_conf -o conf.gro &> /dev/null
gmx_natom=`grep SOL topol.top | awk '{print $2}'`
gmx_natom=`echo "$gmx_natom * $copy_x * $copy_y * $copy_z" | bc`
sed -i "s/SOL.*[0-9]*/SOL $gmx_natom/g" topol.top
$gmx_grompp &> /dev/null
cd ..

# mk jobs script
if [ ! -f $batch_sub ]; then
    echo "cannot find `pwd`/$batch_sub, do nothing, exit"
    exit 1
fi
$tools_dir/mk.batch.sub.sh

cd $cwd
cat $log
mv $log $out_dir

