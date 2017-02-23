#!/bin/bash

targets=$*

for ii in $targets;
do 
    echo doing $ii
    cd $ii
    /home/wanghan/study/ice.melting/templates/tools/cv/locvol --mol-value  -a 	-f out.xtc
    /home/wanghan/study/ice.melting/templates/tools/scripts/mole.smth.traj.sh 	locvol locvol.avg
    /home/wanghan/study/ice.melting/templates/tools/scripts/mole.to.time.py 	locvol.avg/
    /home/wanghan/study/ice.melting/templates/tools/scripts/cv.defect.sh	-f gmx.tool/pbc.xtc -q locvol.avg -t1 .45 -t2 .445 -ot figs.locvol
#    /home/wanghan/study/ice.melting/templates/tools/scripts/plot.script.sh	locvol.avg.def.tcl figs.locvol/
    cd ..
done
# /home/wanghan/study/ice.melting/templates/tools/scripts/plot.hbond.defect.sh
# /home/wanghan/study/ice.melting/templates/tools/scripts/plot.q.defect.sh
