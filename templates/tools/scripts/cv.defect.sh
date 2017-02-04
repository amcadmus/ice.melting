#!/bin/bash

function usage () {
    cat <<EOF
# Detect the defects by the local value of a CV. Consider 2 thresholds to distinguish the metastable state. Generate the plotting script for the VMD.
# optional arguments:
	-h|--help				print this message
	-c|--config file (conf.gro)		tell the gromacs config file. 
	-f|--traj   file (traj.xtc)		tell the gromacs traj file. 
	-q|--loc-q  dir  (steinhardt)		tell the dir of the CV. 
       -t1|--thres1 val  (0.358)		tell the threshold 1 of the CV. 
       -t1|--thres2 val  (0.267)		tell the threshold 2 of the CV. 
       -ot|--o-tcl  dir  (fig.qdef)		output dir for frame tcl script
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

function write_tcl_head () {    
    _my_file_name=$1
    cat > $_my_file_name <<EOF
set mol_id [ mol new $file_conf ]
mol addfile $file_traj waitfor all \$mol_id
set sel [atomselect top all] 
\$sel set radius 0.6
mol addrep \$mol_id
mol addrep \$mol_id
mol modselect 0 \$mol_id not name MW 
mol modselect 1 \$mol_id none
mol modselect 2 \$mol_id none
mol modstyle 0 \$mol_id hbonds
mol modstyle 1 \$mol_id vdw
mol modstyle 2 \$mol_id vdw
mol modcolor 1 \$mol_id ResType
display projection orthographic
display height 4
rotate x by $rot_angle
animate goto 0
EOF
}

function write_tcl_light_head () {    
    _my_file_name=$1
    cat > $_my_file_name <<EOF
EOF
}

file_conf=conf.gro
file_traj=traj.xtc
file_q=steinhardt
thres1=0.358
thres2=0.267
rot_angle=0
out_tcl_dir=figs.qdef

# if [ $# -eq 0 ]; then
#     usage
#     exit
# fi
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
	-t1|--thres1)
	    thres1=$2
	    shift
	    ;;
	-t2|--thres2)
	    thres2=$2
	    shift
	    ;;
	-ot|--o-tcl)
	    out_tcl_dir=$2
	    shift
	    ;;
	-h|--help|*)
	    usage
	    exit
	    ;;
    esac
    shift
done

out_vmd=$file_q.def.tcl
out_file=$file_q.def.out
echo conf: $file_conf
echo traj: $file_traj
echo loc-q: $file_q
echo thres1: $thres1
echo thres2: $thres2
echo out vmd script: $out_vmd
echo out def count : $out_file
echo out frame vmd script dir: $out_tcl_dir

check_file $file_conf
check_file $file_traj
check_folder $file_q
check_folder $out_tcl_dir

targets=`ls $file_q | grep ^step_`
script=$out_vmd

write_tcl_head $script

count=1
rgbs=""
rm -f $out_file
prog1=awkprog1.$$
prog2=awkprog2.$$
echo "{if (\$2 < $thres1 && \$2 > $thres2) {print \$1}}" > $prog1
echo "{if (\$2 < $thres2) {print \$1}}" > $prog2
for ii in $targets; 
do
    output="frame_`printf %09d $count`_$ii"
    step_script="$out_tcl_dir/frame_`printf %09d $count`.tcl"
    if test $count -eq 1; then
	write_tcl_head $step_script
    else
	write_tcl_light_head $step_script
    fi
    echo "animate goto $count" >> $script
    echo "animate goto $count" >> $step_script
#    echo "wait 1" >> $script
    count_frame1=0
    list1=""
    for jj in `grep -v \# $file_q/$ii | awk -f $prog1`
    do
	pjj=$jj
#	pjj=`printf %f $jj | cut -d '.' -f 1`
	if test $count_frame1 -eq 0; then
	    list1="residue $pjj"
	else
	    list1="$list1 or residue $pjj"
	fi
	count_frame1=$(($count_frame1+1))
    done
    count_frame2=0
    list2=""
    for jj in `grep -v \# $file_q/$ii | awk -f $prog2`
    do
	pjj=$jj
#	pjj=`printf %f $jj | cut -d '.' -f 1`
	if test $count_frame2 -eq 0; then
	    list2="residue $pjj"
	else
	    list2="$list2 or residue $pjj"
	fi
	count_frame2=$(($count_frame2+1))
    done
    idx_step=`echo $ii | cut -d '_' -f 2 | cut -d '_' -f 1`
    printf "$ii count $count_frame1 $count_frame2 \r"
    echo "$idx_step $count_frame1 $count_frame2 " >> $out_file
    if test $count_frame1 -ne 0; then
	echo "mol modselect 1 \$mol_id not name MW and ( $list1 )" >> $script
	echo "mol modselect 1 \$mol_id not name MW and ( $list1 )" >> $step_script
    else
	echo "mol modselect 1 \$mol_id none" >> $script
	echo "mol modselect 1 \$mol_id none" >> $step_script
    fi
    if test $count_frame2 -ne 0; then
	echo "mol modselect 2 \$mol_id not name MW and ( $list2 )" >> $script
	echo "mol modselect 2 \$mol_id not name MW and ( $list2 )" >> $step_script
    else
	echo "mol modselect 2 \$mol_id none" >> $script
	echo "mol modselect 2 \$mol_id none" >> $step_script
    fi
    echo "render Tachyon OUT_FOLDER/$output.tachyon" >> $script
    count=$(($count+1))
done
echo ""

rm -fr $prog1 $prog2

# vmd < tmp0.tcl

# for ii in $rgbs;
# do
#     echo convert $ii to png
#     png_name=`echo $ii | sed -e 's/tachyon/png/g'`
# #    convert $ii $png_name
#     tachyon -aasamples 12 $ii -format PNG -o $png_name
# done

# rm -f $rgbs

# cd $out_folder
# count=0; 
# for ii in frame_*png; 
# do 
#     ln -s $ii tmp_`printf %06d $count`.png; 
#     count=$(($count+1)); 
# done
# ffmpeg -i tmp_%06d.png movie.mpg
# cd ..
