GMX_DIR=/home/mi/wanghan/local/gromacs/4.6.7.plumed

source $GMX_DIR/bin/GMXRC

exec_genconf='genconf'
exec_grompp='grompp_mpi'
exec_mdrun="mpirun -n 2 mdrun_mpi"

