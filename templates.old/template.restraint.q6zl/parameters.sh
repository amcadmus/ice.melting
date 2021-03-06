#!/bin/bash

# system setting
system_box="1 1 1"
system_seed_dir=md.seed
system_init_mode=traj
system_init_conf=conf.sol.gro	# in $system_seed_dir
system_init_trr=/home/mi/wanghan/study/ice.melting/templates/template.restraint/traj.trr
system_init_trr_time=27637.821313

# warmup control
system_warming_dir=warming.0.090.0.060
gmx_warm_time=2
gmx_warm_dt=0.001		# in ps
gmx_warm_vel_seed=`date +%s`

# run time control
system_running_dir=running.0.090.0.060
gmx_time=2000			# in ps
gmx_dt=0.002			# in ps
gmx_temperature=252.10		# in K
gmx_tau_t=0.1			# in ps
gmx_tau_box=2			# in ps
gmx_pcoupltype=isotropic
plumed_q6_at=0.060000
plumed_q6_kappa=1000000		# take 432 water system as reference
plumed_zl_at=0.060000
plumed_zl_kappa=1000000		# take 432 water system as reference
plumed_zl_ref_numb=6
plumed_zl_comp_shift=true
plumed_nlist_rc=0.7
plumed_nlist_feq=50

# output control
gmx_xtc_feq=1.0			# in ps
gmx_trr_feq=0			# in ps
gmx_energy_feq=1.0		# in ps

# commands
grompp_cmd=grompp_mpi
mdrun_cmd=mdrun_mpi
#gmx_dir=$HOME/local/gromacs/4.6.5


