#!/bin/bash

temperature=200
nsteps=100000

sed -e "s/ref_t.*=.*/ref_t = $temperature/g" grompp.mdp | \
sed -e "s/nsteps.*=.*/nsteps = $nsteps/g" > tmp.tmp
mv -f tmp.tmp grompp.mdp

grompp_mpi
mpirun -n 4 mdrun_mpi  -v 
echo 3 | trjconv_mpi -f traj.xtc -o oxygen.xtc -n index.ndx -b 50
mpirun -n 4 plumed driver --mf_xtc oxygen.xtc --plumed plumed.input

