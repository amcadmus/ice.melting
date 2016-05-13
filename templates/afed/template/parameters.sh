# system setting
system_seed_dir=md
system_start_conf=conf.sol.gro	# in $system_seed_dir
system_conf_copies=2,2,2

# run time control
md_time=100000
md_dt=0.002
md_rand_seed=`date +%s`
md_T=252.1
md_tau_t=0.2
md_tau_p=1.0
md_ener_freq=0.2
md_traj_freq=2.0
afed_T=2.0e4
afed_tau=10

gmx_bin_dir=$HOME/local/gromacs/4.6.5/bin
gmx_genconf=genconf
gmx_grompp=grompp

