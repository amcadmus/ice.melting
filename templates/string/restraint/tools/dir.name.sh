#!/bin/bash

source parameters.sh
source env.sh

cwd=`pwd`

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
system_template_conf=$cwd/$system_seed_dir/conf.${system_start_state}.gro
system_template_copies=$system_conf_copies

# copy 
copy_x=`echo $system_template_copies | cut -d ',' -f 1`
copy_y=`echo $system_template_copies | cut -d ',' -f 2`
copy_z=`echo $system_template_copies | cut -d ',' -f 3`
copy_x=`printf %02d $copy_x`
copy_y=`printf %02d $copy_y`
copy_z=`printf %02d $copy_z`
natom=`wc -l $system_template_conf | awk '{print $1}' `
natom=`echo "($natom - 3) / 4 * $copy_x * $copy_y * $copy_z" | bc`
natom=`printf %07d $natom`

# parse constraints.
top_file=tops/top.input.$system_top_style
numb_res=`grep HarmonicRestraint $top_file | wc -l`
numb_centers=`echo $res_centers | awk -F ',' '{print NF}'`
res_name=""
for ii in `seq 1 $numb_centers`
do
    key=`grep HarmonicRestraint $top_file | head -n $ii | tail -n 1 | awk '{print $1}'`
    center_value=`echo $res_centers | cut -d ',' -f $ii`
    center_value=`printf %.04f $center_value`
    res_name=$res_name.${key}$'-'${center_value}
done

if [ ${#dir_prefix} -gt 0 ]; then
    my_dir_prefix=${dir_prefix}'.'
else
    my_dir_prefix=${dir_prefix}
fi
out_dir=${my_dir_prefix}${system_start_state}.${natom}$res_name

echo $out_dir

