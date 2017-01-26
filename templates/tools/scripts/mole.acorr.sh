#!/bin/bash

moasp_dir=~/study/moasp/build/src/app/

if test $# -lt 1 ; then
    echo "# compute the auto-correlation function for mole traj files matches target/mol_*"
    echo "# usage"
    echo $0 target [upper_time_bound]
    exit
fi
upper=1000
if test $# -ge 2; then
    upper=$2
fi

folder=$1

echo "# do folder $folder"
echo "# acc time up to $upper"

cd $folder 
rm -f acorr.all.out
for ii in mol_*; 
do 
    printf "# doing %s \r" $ii
    result=`$moasp_dir/msp_acorr -f $ii -t 2 -d 3 -u $upper --scale -o acorr.$ii | grep -v \# | grep -v ^$`; 
    aidx=`echo $ii | cut -d '_' -f 2`;
    echo $aidx $result >> acorr.all.out ; 
done
echo ""
