#!/bin/bash

if test $# -lt 1; then
    echo usage
    echo $0 [ string [ string [ string... ] ] ]
    exit
fi

targets=$*

for ii in $targets
do 
    if test ! -d $ii; then
	echo "no string $ii"
	continue
    fi
    nodes=`ls $ii | grep node.0`
    for jj in $nodes
    do
	node_dir=$ii/$jj
	if test -f $node_dir/out.png; then
	    echo existing $node_dir/out.png, not plotting
	    continue
	fi
	if test ! -f $node_dir/out.gro; then
	    echo "no gro file $node_dir/out.gro"
	    exit
	fi
	cd $node_dir	
	tmptcl=tmptcl.$$
	cat > $tmptcl <<EOF
mol new out.gro 
display projection orthographic
display height 4
rotate x by 0
render Tachyon tmp.tachyon
EOF
	vmd < $tmptcl
	sleep 1
	tachyon -aasamples 0 tmp.tachyon -format PNG -o out.png
#	rm -f tmp.tachyon $tmptcl
	cd ../..
    done
done
	
	      
