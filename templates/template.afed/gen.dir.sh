#!/bin/bash

source parameters.sh

dir_base=`pwd`

dir_current=$system_running_dir
if test -d $dir_current; then
    echo "exist dir $dir_current, continue"
    exit
fi
if test ! -d md.seed; then
    echo "the dir md.seed does not exist, exit"
    exit
fi

echo "# generate dir $dir_current"
cp -a md.seed $dir_current    

cd $dir_current

echo "# gen conf"
genconf -nbox $system_box -f $system_base_conf -o out.gro
mv -f out.gro conf.gro

echo "# update top"
nline=`wc -l conf.gro | awk '{print $1}'`
nline=$(($nline-3))
nmol=`echo "$nline / 4" | bc`
echo "## the system is estimated to have $nmol tip4p water molecules"
sed "s/SOL.*/SOL $nmol/g" topol.top > tmp.top
mv -f tmp.top topol.top

cd $dir_base
./update.dir.sh

