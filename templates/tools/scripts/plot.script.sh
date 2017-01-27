#!/bin/bash

if test $# -ne 2; then
    echo usage
    echo $0 tcl_script out_folder
    exit
fi

tcl_script=$1
out_folder=$2
if test ! -d $out_folder; then
    mkdir -p $out_folder
fi

sed -e "s/OUT_FOLDER/$out_folder/g" $tcl_script > tmp.$tcl_script

vmd < tmp.$tcl_script
tachyon_files=`grep $out_folder tmp.$tcl_script | grep tachyon | awk '{print $3}'`
rm -f tmp.$tcl_script

for ii in $tachyon_files;
do
    echo convert $ii to png
    png_name=`echo $ii | sed -e 's/tachyon/png/g'`
#    convert $ii $png_name
    tachyon -aasamples 12 $ii -format PNG -o $png_name
done

cd $out_folder
count=0; 
for ii in frame_*png; 
do 
    ln -s $ii tmp_`printf %06d $count`.png; 
    count=$(($count+1)); 
done
ffmpeg -i tmp_%06d.png movie.mpg
cd ..


