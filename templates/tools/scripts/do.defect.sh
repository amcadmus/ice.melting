#!/bin/bash

targets=$*

for ii in $targets;
do 
    echo doing $ii
    cd $ii
    if test ! -f gmx.tool/pbc.xtc; then
	cd gmx.tool
	ln -sf ../conf.gro .
	grompp 
	echo 0 | trjconv -f ../out.xtc -pbc mol  -o pbc.xtc
	cd ..
    fi
#    /home/wanghan/study/ice.melting/templates/tools/h.bond/hbond  --detail --mol-defect		-f gmx.tool/pbc.xtc -t 2
    /home/wanghan/study/ice.melting/templates/tools/cv/steinhardt --detail --mol-value  -a 	-f gmx.tool/pbc.xtc -t 2
    /home/wanghan/study/ice.melting/templates/tools/cv/locvol     --detail --mol-value  -a 	-f gmx.tool/pbc.xtc -t 2
#    /home/wanghan/study/ice.melting/templates/tools/scripts/mole.acorr.sh 	nhbond 1000
    /home/wanghan/study/ice.melting/templates/tools/scripts/mole.smth.traj.sh 	steinhardt steinhardt.avg 
    /home/wanghan/study/ice.melting/templates/tools/scripts/mole.to.time.py 	steinhardt.avg 
    /home/wanghan/study/ice.melting/templates/tools/scripts/mole.smth.traj.sh 	locvol locvol.avg 
    /home/wanghan/study/ice.melting/templates/tools/scripts/mole.to.time.py 	locvol.avg 
    /home/wanghan/study/ice.melting/templates/tools/scripts/cv.defect.sh 	-f gmx.tool/pbc.xtc -q steinhardt.avg -t1 .324 -t2 .274
#    /home/wanghan/study/ice.melting/templates/tools/scripts/cv.crystal.sh 	-f gmx.tool/pbc.xtc -q steinhardt.avg -t1 .324 -t2 .274
#    /home/wanghan/study/ice.melting/templates/tools/scripts/plot.script.sh	steinhardt.avg.def.tcl figs.qdef/
    cd ..
done
# /home/wanghan/study/ice.melting/templates/tools/scripts/plot.hbond.defect.sh
# /home/wanghan/study/ice.melting/templates/tools/scripts/plot.q.defect.sh
