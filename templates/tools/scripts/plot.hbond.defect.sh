#!/bin/bash

function usage () {
    cat <<EOF
# plot the defects defined by H-bonds via VMD
# optional arguments:
	-h|--help			print this message
	-c|--config file (conf.gro)	tell the gromacs config file. 
	-f|--traj   file (traj.xtc)	tell the gromacs traj file. 
	-b|--hbond  dir  (nhbond)	tell the dir of h-bond info. 
	-o|--output dir  (figs.hdef)	the output dir.	
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
file_hbond=nhbond
out_folder=figs.hdef
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
	-b|--hbond)
	    file_hbond=$2
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
        # *)
        #     echo "Unknow option $1, the program dies"
        #     exit
    esac
    shift
done

echo conf: $file_conf
echo traj: $file_traj
echo hbond: $file_hbond
echo out: $out_folder

check_file $file_conf
check_file $file_traj
check_folder $file_hbond
mkdir -p $out_folder

targets=`ls $file_hbond | grep time_`
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
for ii in $targets; 
do
    output="frame_`printf %09d $count`_$ii"
    echo "animate goto $count" >> $script
#    echo "wait 1" >> $script
    count_frame=0
    list=""
    # for jj in `grep -v '2 2' $file_hbond/$ii | grep -v '2 3' | grep -v '3 2' | grep -v \# | awk '{print $1}'` 
    for jj in `grep -v '2 2' $file_hbond/$ii | grep -v \# | awk '{print $1}'` 
    do
	if test $count_frame -eq 0; then
	    list="residue $jj"
	else
	    list="$list or residue $jj"
	fi
	count_frame=$(($count_frame+1))
    done
    echo "mol modselect 1 \$mol_id $list" >> $script
#    echo "wait 0.5" >> $script
#    echo "render snapshot $out_folder/$output.rgb" >> $script
    echo "render Tachyon $out_folder/$output.tachyon" >> $script
    rgbs="$rgbs $out_folder/$output.tachyon"
    count=$(($count+1))
done

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

# rot_angle=0
# if test $# -lt 2 ; then
#     echo usage
#     echo "$0 <input_conf> <output_fig (no extension)> [rot_angle_by_x]"
#     exit
# fi
# if test $# -ge 3 ; then
#     rot_angle=$3
# fi

# input=$1
# output=$2

# cat > tmp0.tcl <<EOF
# mol new $input
# display projection orthographic
# rotate x by $rot_angle
# render snapshot $output.rgb
# EOF
# vmd < tmp0.tcl

# convert $output.rgb $output.png

