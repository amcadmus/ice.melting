#!/usr/bin/env python3

import os
import argparse
import re
import numpy as np
import logging
import subprocess as sp
import enum as Enum
from StringForce import StringForce
from subprocess import Popen, PIPE
import tube_genjobs as tbgen

import glob
from BatchJob import JobStatus
from PBSJob import PBSJob
from SlurmJob import SlurmJob
import time
import datetime

def submit_step (step_name,
                 max_numb_job) :    
    base_dir = os.getcwd()
    sf = StringForce ("template.string")
    if False == os.path.exists (step_name) :
        raise RuntimeError ("no step dir " + step_name + " should be generated before submitting")
    string = np.loadtxt (step_name + "/string.out")
    os.chdir (step_name)
    step_dir = os.getcwd() + "/"
    job_list = []
    fpenv = open("env.sh", "r")
    for line in fpenv :
        line = line.rstrip()
        if re.search ("batch_system=", line) :
            batch_system = line.split("=")[-1]
            break
    for node_idx in range (string.shape[0]) :            
        node_name = sf.mk_node_name (node_idx)
        os.chdir (node_name)
        node_dir = os.getcwd() + "/"
        sp.check_call ([step_dir + sf.cmd_job_scpt])
        if batch_system == "PBS" :
            job = PBSJob (node_dir, "PBS.sub")
        elif batch_system == "Slurm" :
            job = SlurmJob (node_dir, "Slurm.sub")
        else :
            raise RuntimeError ("Unknown batch system")
        job_id = job.submit ()
        if job.check_status () != JobStatus.finished :
            job_list.append (job)
        os.chdir (step_dir)
        if len(job_list) >= max_numb_job :
            break
    os.chdir (base_dir)        
    return job_list

def wait_step (job_list) :
    while True :
        find_unfinish = False
        for job in job_list :
            stat = job.check_status ()
            if stat == JobStatus.terminated :
                logging.info ("find terminated job %s, wait and check again " % job.get_job_id())
                time.sleep (10)
                stat = job.check_status ()
                if stat == JobStatus.terminated :
                    msg = "find terminated job %s. exit. should restart" % job.get_job_id()
                    logging.error (msg)
                    raise RuntimeError (msg)
                else :
                    logging.info ("check passed")
            if stat != JobStatus.finished :
                find_unfinish = True
                break
        if find_unfinish == False :
            return
        else :
            time.sleep (60)
        print ("# checked at " + time.strftime("%Y-%m-%d %H:%M:%S"))            
    
    
def exec_jobs (max_numb_job) :
    dirnames = glob.glob ("step.*")
    dirnames.sort()
    if os.path.exists ("tag_fin_step") :
        fin_steps = [ line.rstrip('\n') for line in open("tag_fin_step") ]
        step_list = [ step for step in dirnames if step not in fin_steps ]
    else :
        fp = open ("tag_fin_step", "w")
        fp.write (tbgen.make_dir_name(0) + "\n")
        fp.close ()
        step_list = [ step for step in dirnames if not step == tbgen.make_dir_name(0) ]

    logging.info ("will do steps %s" % step_list)

    for ii in step_list :
        while True :
            job_list = submit_step (ii, max_numb_job)
            if len(job_list) == 0 :
                break
            wait_step (job_list)            
        fp = open ("tag_fin_step", "a")
        fp.write (ii + "\n")
        fp.close ()
        logging.info ("finished step %s" % ii)

def main ():
    parser = argparse.ArgumentParser(
        description="*** Execute the jobs of a tube. ***")
    parser.add_argument('-m', '--max-numb-jobs', type=int, default = 100,
                        help='The maximum number of allowed jobs. ')
    args = parser.parse_args()

    logging.basicConfig (filename="tube_exec.log", filemode="w", level=logging.INFO, format='%(asctime)s %(message)s')

    exec_jobs (args.max_numb_jobs)        

if __name__ == "__main__":
    main ()
        
