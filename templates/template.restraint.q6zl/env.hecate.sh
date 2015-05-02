GMX_DIR=/home/han.wang/local/gromacs/4.6.7.plumed

source $GMX_DIR/bin/GMXRC

exec_genconf='srun -n 1 genconf'
exec_grompp='srun -n 1 grompp_mpi'
exec_mdrun="srun -n 6 mdrun_mpi"
