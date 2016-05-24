#!/usr/bin/env python3

import os
import logging
import sys
from enum               import Enum
from subprocess         import Popen, PIPE
Logger = logging.getLogger("string_method")

class JobStatus (Enum) :
    unsubmitted = 1
    waiting = 2
    running = 3
    terminated = 4
    finished = 5
    unknow = 100

class BatchJob (object):
    """
    Abstract class of a batch job
    It submit a job (leave the id in file tag_jobid)
    It check the status of the job (return JobStatus)
    NOTICE: I assume that when a job finishes, a tag file named tag_finished should be touched by the user.
    TYPICAL USAGE:
    job = DERIVED_BatchJob (dir, script)
    job.submit ()
    stat = job.check_status ()    
    """
    def __init__ (self,                         
                  job_dir = "",                         # dir of the job
                  job_script = "",                      # name of the job script
                  job_finish_tag = "tag_finished",      # name of the tag for finished job
                  job_id_file = "tag_jobid") :          # job id if making an existing job
        self.job_dir = job_dir
        self.job_script = job_script
        self.job_id_file = job_dir + "/" + job_id_file
        self.job_finish_tag = job_dir + "/" + job_finish_tag
        self.cwd = os.getcwd()
        self.submit_cmd = str(self.submit_command())
    def get_job_id (self) :
        if True == os.path.exists (self.job_id_file) :
            fp = open (self.job_id_file, 'r')
            job_id = fp.read ()
            return str(job_id)
        else :
            return ""
    def submit_command (self) :
        """
        submission is 
        $ [command] [script]
        """
        Logger.error("submit_command not implemented")
    def check_status (self):
        Logger.error ("check_status not implemented")
    def submit (self) :
        if self.get_job_id () != "" :
            stat = self.check_status()
            if stat != JobStatus.terminated :
                if stat == JobStatus.unknow :
                    raise RuntimeError ("unknown job status, terminate!")
                Logger.warning ("# job " +
                                self.get_job_id() +
                                " already submitted (waiting, running or finished), would not submit again")
                return self.get_job_id()
            else :
                Logger.warning ("# find terminated job " + self.get_job_id() + ", submit again")                
        if (False == os.path.isdir (self.job_dir) ) :
            Logger.error ("cannot find job dir " + self.job_dir)
        abs_job_script = self.job_dir + "/" + self.job_script
        if False == os.path.exists (abs_job_script) :
            Logger.error ("cannot find job script " + abs_job_script)
        cwd = os.getcwd()
        os.chdir (self.job_dir)
        ret = Popen([self.submit_cmd, self.job_script], stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if str(stderr, encoding='ascii') != "":
            Logger.error (stderr)
        job_id = str(stdout, encoding='ascii').replace('\n','')
        print ("# submitted job has ID: " + job_id)
        fp = open (self.job_id_file, 'w')
        fp.write (job_id)
        fp.close()
        os.chdir (cwd)
        return self.get_job_id()
        
class PBSJob (BatchJob) :
    def submit_command (self):
        return "qsub"
    def check_status (self):
        job_id = self.get_job_id ()
        if len(job_id) == 0 :
            return JobStatus.unsubmitted
        ret = Popen (["qstat", job_id],  stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if (ret.returncode != 0) :
            if str("qstat: Unknown Job Id") in str(stderr, encoding='ascii') :
                if os.path.exists (self.job_finish_tag) :
                    return JobStatus.finished
                else :
                    return JobStatus.terminated
            else :
                Logger.error ("status command " + "qstat" + " fails to execute")
                Logger.error ("erro info: " + str(stderr, encoding='ascii'))
                Logger.error ("return code: " + str(ret.returncode))
                sys.exit ()
        status_line = str(stdout, encoding='ascii').split ('\n')[-2]
        status_word = status_line.split ()[-2]        
#        print (status_word)
        if      status_word in ["Q","H"] :
            return JobStatus.waiting
        elif    status_word in ["R"] :
            return JobStatus.running
        elif    status_word in ["C","E","K"] :
            if os.path.exists (self.job_finish_tag) :
                return JobStatus.finished
            else :
                return JobStatus.terminated
        else :
            return JobStatus.unknown        

# if __name__ == "__main__" :
#     job = PBSJob ("./test", "cu01.sleep")
#     job.submit ()
#     print ("submit done")
#     stat = job.check_status ()
#     print ("check done")
#     print (stat)
    
