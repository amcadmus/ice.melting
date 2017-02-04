#!/bin/bash 

tq6=0.324
tdispl=0.4

echo "# usage "
echo "# $0 [threshold_q6] [threshold_displ]"

if test $# -ge 1; then
    tq6=$1
fi

if test $# -ge 2; then
    tdispl=$2
fi

echo "# threshold_q6:	 $tq6"
echo "# threshold_dipl:	 $tdispl"

echo "# index q6 displ"

for ii in `seq 0 3455`
do 
    pii=`printf %06d $ii`
    q6=`tail -n 1 steinhardt/mol_$pii | awk '{print $3}'`
    displ=`tail -n 1 displ/mol_$pii | awk '{print $3}'`
#    echo $ii $q6 $displ
    if [ `echo "$q6 > $tq6" | bc -l` -eq 1 ] && [ `echo "$displ > $tdispl" | bc -l` -eq 1 ]; then
	echo $ii $q6 $displ
    fi
done

