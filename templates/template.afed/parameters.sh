#!/bin/bash

# system setting
system_box="1 1 1"
system_seed_dir=md.seed
system_base_conf=conf.gro	# in $system_seed_dir
system_running_dir=running

# run time control
gmx_time=1000000		# in ps
gmx_dt=0.002			# in ps
gmx_temperature=252.10		# in K
gmx_tau_t=0.1			# in ps
gmx_tau_box=40			# in ps
gmx_pcoupltype=anisotropic
plumed_cv_temperature=100000	# in K
plumed_cv_tau=1000		# in ps
plumed_cv_kappa=500000		# take 432 water system as reference
plumed_nlist_rc=0.7
plumed_nlist_feq=50

# output control
gmx_xtc_feq=10.0		# in ps
gmx_trr_feq=0			# in ps
gmx_energy_feq=10		# in ps

# commands
grompp_cmd=grompp_mpi
mdrun_cmd=mdrun_mpi
#gmx_dir=$HOME/local/gromacs/4.6.5


