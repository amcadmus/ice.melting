#!/bin/bash

source parameters.sh
source env.sh
source $GMX_DIR/bin/GMXRC

dir_base=`pwd`
gmx_nsteps=`echo "($gmx_dt*0.5+$gmx_equi_time) / $gmx_dt" | bc`

nproc=`echo $temperatures | wc -w`

for temperature in $temperatures ;
do
    dir_current=dir.$temperature
    if test ! -d $dir_current; then
	echo "does not exist dir $dir_current, exit"
	exit
    fi
    echo "# grompp in dir $dir_current"
    cd $dir_current
    $grompp_cmd -n index.ndx &> grompp.log
    cd $dir_base
done

command="mpirun -n $nproc $mdrun_cmd -multidir dir* -v -nsteps $gmx_nsteps"
echo "#run with command: $command"
$command

for temperature in $temperatures ;
do
    dir_current=dir.$temperature
    if test ! -d $dir_current; then
	echo "does not exist dir $dir_current, exit"
	exit
    fi
    echo "# update config $dir_current"
    cd $dir_current
    mv -f confout.gro conf.gro
    cd $dir_base
done

# for temperature in $temperatures ;
# do 
#     dir_current=dir.$temperature
#     if test ! -d $dir_current; then
# 	echo "does not exist dir $dir_current, exit"
# 	exit
#     fi    
#     echo "# run in dir $dir_current"
#     cd $dir_current
#     $grompp_cmd -n index.ndx
#     $mdrun_equi_cmd -nsteps $gmx_nsteps
#     mv -f confout.gro conf.gro
#     cd $dir_base
# done
    
