#!/bin/bash

source env.sh
source parameters.sh
tools_dir=`cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd`

dep_job_id=""
if test $# -ge 1; then
    dep_job_id=$1
fi

nnode=`echo "$numb_proc / $numb_proc_per_node" | bc`
rest_proc=`echo "$numb_proc - $nnode * $numb_proc_per_node" | bc `
if [ $rest_proc -gt 0 ]; then
    if [ $nnode -eq 0 ]; then
	nnode=$(($nnode+1))
	numb_proc_per_node=$rest_proc
    else
	echo "cannot make pbs job with $numb_proc for $nnode nodes and $numb_proc_per_node proc per node"
    fi
fi    

echo "#PBS -S /bin/bash"			>  $batch_sub
echo "#PBS -N moasp"				>> $batch_sub
echo "#PBS -q $batch_queue"			>> $batch_sub
echo "#PBS -l nodes=$nnode:ppn=$numb_proc_per_node"	 >> $batch_sub
echo "#PBS -l walltime=$job_hour:$job_min:00"	>> $batch_sub
if [ ${#dep_job_id} -ne 0 ]; then
    echo "#PBS -W depend=afterok:$dep_job_id"	>> $batch_sub
fi
echo "#PBS -l cput=$cput_hour:00:00"		>> $batch_sub
echo "#PBS -j oe"				>> $batch_sub
echo ""						>> $batch_sub
echo 'NP=`cat $PBS_NODEFILE | wc -l`'		>> $batch_sub
echo 'source /opt/intel/composer_xe_2015/bin/compilervars.sh intel64'	>> $batch_sub
echo 'export JASMIN_NUM_SMP_THREADS=1'		>> $batch_sub
echo 'export JASMIN_NUM_DSM_THREADS=1'		>> $batch_sub
echo 'export OMP_NUM_THREADS=1'			>> $batch_sub
echo 'echo np is $NP'				>> $batch_sub
echo ""						>> $batch_sub
echo 'cd $PBS_O_WORKDIR'			>> $batch_sub
echo "$moasp_bin_dir/$moasp_ipp -n \$NP"	>> $batch_sub
echo 'if [ $? -ne 0 ]; then exit 1; fi'		>> $batch_sub
echo ""						>> $batch_sub
echo "$mpirun_command -genv I_MPI_DEVICE rdma -hostfile \$PBS_NODEFILE -n \$NP $moasp_bin_dir/$moasp_run -v &> log" >> $batch_sub
echo 'if [ $? -ne 0 ]; then exit 1; fi'		>> $batch_sub
echo ""						>> $batch_sub
echo "$tools_dir/last.conf.sh"			>> $batch_sub
echo 'if [ $? -eq 0 ]; then touch tag_finished; fi' >> $batch_sub
echo ""						>> $batch_sub

