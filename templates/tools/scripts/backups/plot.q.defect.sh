#!/bin/bash

function usage () {
    cat <<EOF
# plot the defects defined by Steinhardt parameter via VMD
# optional arguments:
	-h|--help			print this message
	-c|--config file (conf.gro)	tell the gromacs config file. 
	-f|--traj   file (traj.xtc)	tell the gromacs traj file. 
	-q|--loc-q  dir  (steinhardt)	tell the dir of local Steinhardt parameter. 
	-t|--thres  val  (0.358)	tell the threshold of Steinhardt parameter. 
	-o|--output dir  (figs.qdef)	the output dir.	
EOF
}

function check_file () {
    file=$1
    if test ! -f $file; then
	echo no file $file
    fi
}

function check_folder () {
    folder=$1
    if test ! -d $folder; then
	echo no file $folder
    fi
}

file_conf=conf.gro
file_traj=traj.xtc
file_q=steinhardt
thres=0.358
out_file=$file_q.def.out
out_folder=figs.qdef
rot_angle=0

if [ $# -eq 0 ]; then
    usage
    exit
fi
while [ $# -ge 1 ]; 
do
    key=$1
    case $key in
	-c|--config)
	    file_conf=$2
	    shift
	    ;;
	-f|--traj)
	    file_traj=$2	    
	    shift
	    ;;
	-q|--loc-q)
	    file_q=$2
	    shift
	    ;;
	-t|--thres)
	    thres=$2
	    shift
	    ;;
	-o|--output)
	    out_folder=$2	    
	    shift
	    ;;
	-h|--help|*)
	    usage
	    exit
	    ;;
    esac
    shift
done

echo conf: $file_conf
echo traj: $file_traj
echo loc-q: $file_q
echo thres: $thres
echo out: $out_folder
echo out file: $out_file

check_file $file_conf
check_file $file_traj
check_folder $file_q
mkdir -p $out_folder

targets=`ls $file_q | grep ^step_`
script=tmp0.tcl

cat > $script <<EOF
set mol_id [ mol new $file_conf ]
mol addfile $file_traj waitfor all \$mol_id
mol addrep \$mol_id
mol modselect 0 \$mol_id all
mol modselect 1 \$mol_id none
mol modstyle 0 \$mol_id hbonds
mol modstyle 1 \$mol_id vdw
display projection orthographic
rotate x by $rot_angle
animate goto 0
EOF

count=1
rgbs=""
rm -f $out_file
echo "{if (\$2 < $thres) {print \$1}}" > awkprog
for ii in $targets; 
do
    output="frame_`printf %09d $count`_$ii"
    echo "animate goto $count" >> $script
#    echo "wait 1" >> $script
    count_frame=0
    list=""
    for jj in `grep -v \# $file_q/$ii | awk -f awkprog`
    do
	if test $count_frame -eq 0; then
	    list="residue $jj"
	else
	    list="$list or residue $jj"
	fi
	count_frame=$(($count_frame+1))
    done
    idx_step=`echo $ii | cut -d '_' -f 2 | cut -d '_' -f 1`
    printf "$ii count $count_frame \r"
    echo "$idx_step $count_frame " >> $out_file
    if test $count_frame -ne 0; then
	echo "mol modselect 1 \$mol_id $list" >> $script
    fi
#    echo "wait 0.5" >> $script
#    echo "render snapshot $out_folder/$output.rgb" >> $script
#    echo "render PostScript $out_folder/$output.eps" >> $script
    echo "render Tachyon $out_folder/$output.tachyon" >> $script
    rgbs="$rgbs $out_folder/$output.tachyon"
    count=$(($count+1))
done
echo ""

vmd < tmp0.tcl

for ii in $rgbs;
do
    echo convert $ii to png
    png_name=`echo $ii | sed -e 's/tachyon/png/g'`
#    convert $ii $png_name
    tachyon -aasamples 12 $ii -format PNG -o $png_name
done

rm -f $rgbs

cd $out_folder
count=0; 
for ii in frame_*png; 
do 
    ln -s $ii tmp_`printf %06d $count`.png; 
    count=$(($count+1)); 
done
ffmpeg -i tmp_%06d.png movie.mpg
cd ..
