#!/bin/bash

source parameters.sh
source env.sh
source $gmx_bin_dir/GMXRC.bash

tools_dir=$(cd ${0%/*} && echo $PWD)
cwd=`pwd`
log=gen.dir.log
rm -f $log

if test $# -lt 2 ; then
    echo usage:
    echo $0 out_dir init_folder
    exit 1   
fi
out_dir=$1
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

init_folder=$2
if test ! -d $init_folder; then
    echo "no init folder $init_folder, do nothing"
    exit 1
fi
system_ref_conf=confs/$system_start_conf
system_start_conf=$init_folder/out.gro

traj_n_freq=`echo "($md_traj_freq + 0.5 * $md_dt) / $md_dt " | bc`
ener_n_freq=`echo "($md_ener_freq + 0.5 * $md_dt) / $md_dt " | bc`
natom=`wc -l $system_ref_conf | awk '{print $1}' `
natom=`echo "($natom - 3)" | bc`
nmol=`echo "$natom / 4" | bc`

# parse constraints.
top_file=tops/top.input.$system_top_style
numb_res=`grep "HarmonicRestraint\|VolumeConstraint" $top_file | wc -l`
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
    key=`grep "HarmonicRestraint\|VolumeConstraint" $top_file | head -n $ii | tail -n 1 | awk '{print $1}'`
    center_value=`echo $res_centers | cut -d ',' -f $ii`
    center_value=`printf %.04f $center_value`
    res_name=$res_name.${key}$'-'${center_value}
done

# make output dir
if test -d $out_dir; then
    echo "existing dir $out_dir, do nothing"
    exit 10
fi

# copy dir
echo "# generate dir $out_dir"			>> $log
echo "# using conf file $system_start_conf"	>> $log
echo "# using top file $top_file"		>> $log
cp -a $system_seed_dir $out_dir
cp $top_file $out_dir/top.input
cp parameters.sh $out_dir
cp -L env.sh $out_dir
if test ${#init_folder} -gt 0; then
    echo "$init_folder/" > $out_dir/job_dep
fi

$tools_dir/update.sh $out_dir

cd $out_dir

# prepare conf
ln -s $system_start_conf conf.gro

# update gromacs tools
cd gmx.tool
ln -s ../../$system_ref_conf ./conf.gro
sed -i "s/SOL.*[0-9]*/SOL $nmol/g" topol.top
$gmx_grompp &> /dev/null
cd ..

cd $cwd
cat $log
mv $log $out_dir

