#!/bin/bash

if test $# -ne 1; then
    echo usage
    echo $0 ring_dir
    exit
fi

idir=$1
odir=vmd.$idir

if test ! -d $idir; then
    echo no dir $idir, no nothing
    exit
fi
if test ! -d $odir; then
    mkdir -p $odir
fi

cwd=`pwd`
cd $idir
nframe=`ls def_ring* | wc -l`
cd $cwd

for ii in `seq 1 $nframe`
do
    pi=`printf %09d $ii`
    ofile=$odir/frame.$pi.tcl
    cat > $ofile << EOF    
set mol_id [ mol new conf.gro ]
mol addfile gmx.tool/pbc.xtc waitfor all \$mol_id
set sel [atomselect top all] 
\$sel set radius 0.6
mol modselect 0 \$mol_id not name MW 
mol modstyle 0 \$mol_id hbonds
display projection orthographic
display height 4
rotate x by 0
animate goto 0
animate goto $ii
mol rep hbonds 3.0 40 3
EOF
    ring_file=$idir/`ls $idir | grep def_ring | head -n $ii | tail -n 1`
    nlines=`cat $ring_file | wc -l`
    echo $ring_file $nlines
    for jj in `seq 1 $nlines`
    do
	nwords=`head -n $jj $ring_file | tail -n 1 | wc -w`
#	echo "mol color ColorID $nwords" >> $ofile
	echo "mol addrep \$mol_id" >> $ofile
    done	      
    for jj in `seq 1 $nlines`
    do
	echo "mol modselect $jj \$mol_id none" >> $ofile
    done
    for jj in `seq 1 $nlines`
    do
	nwords=`head -n $jj $ring_file | tail -n 1 | wc -w`
	if test $nwords -eq 5; then
	    echo "mol modstyle $jj \$mol_id vdw" >> $ofile
	else
#	    echo "mol modstyle $jj \$mol_id hbonds" >> $ofile
	    echo "mol modstyle $jj \$mol_id vdw" >> $ofile
	fi
    done	      
    
    for jj in `seq 1 $nlines`
    do
	line=`head -n $jj $ring_file | tail -n 1`
	list=""
	count=0
	for kk in $line
	do
	    if test $count -eq 0; then
		list="residue $kk"
	    else 
		list="$list or residue $kk"
	    fi
	    count=$((count+1))
	done
	echo "mol modselect $jj \$mol_id not name MW and ( $list )" >> $ofile
    done
done

