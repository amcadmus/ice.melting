#!/usr/bin/env python3

import os
import sys
import logging
Logger = logging.getLogger("string_method")
from subprocess import Popen, PIPE
import subprocess as sp
import re
import numpy as np
from PBSJob import PBSJob
from PBSJob import JobStatus
import time
import datetime

class StringForce (object) :
    def __init__ (self,
                  string_template  = "template.string",
                  dep_sec = 4,
                  ) :
        self.base_dir = os.getcwd() + '/'
        self.string_template = string_template
        self.cmd_gen_dir = "tools/gen.dir.sh"
        self.cmd_update_dir = "tools/update.sh"
        self.cmd_job_scpt = "tools/mk.batch.sub.sh"
        self.cmd_stat_scpt = "tools/cmpt.f.sh"
        self.dep_sec_size = dep_sec
        
        if False == os.path.exists (self.string_template) :
            Logger.error ("cannot find string template " + self.string_template)
        os.chdir (self.string_template)
        testfiles = []
        testfiles.append ("parameters.sh")
        testfiles.append (self.cmd_gen_dir)
        testfiles.append (self.cmd_job_scpt)
        for tfile in testfiles :
            if False == os.path.exists (tfile) :
                Logger.error ("no file " +  tfile)
                sys.exit ()
        os.chdir (self.base_dir)
    
    def compute (self,
                 step,
                 string,
                 ) :
        self.generate_string (step, string)
        job = self.submit_string (step)
        self.wait_string (job)
        self.write_tag (step)
        force = self.statistic_string (step)
        return (force)

    def write_tag (self,
                   step,
                   ) :
        fp = open ("tag_fin_string", "a")
        fp.write (str(step) + "\n")
        fp.close ()

    def statistic_string (self,
                          step,
                          ):
        string_name = self.mk_string_name (step)
        if False == os.path.exists (string_name) :
            raise RuntimeError ("no string dir " + string_name + " should be generated before submitting")
        os.chdir (string_name)
        string = np.loadtxt ("string.out")
        sp.check_call (["rm -f force.out"], shell=True)
        for node_idx in range (string.shape[0]) :
            node_name = self.mk_node_name (node_idx)
            sp.check_call ([self.cmd_stat_scpt + " " + node_name + " | grep -v \# >> force.out"], shell = True)
        force = np.loadtxt ("force.out")
        sp.check_call (["paste string.out force.out > all.out"], shell = True)
        os.chdir (self.base_dir)
        return force

    def wait_string (self,
                     job_list,
                     ) :
        while True :
            find_unfinish = False
            for job in job_list :
                stat = job.check_status ()
                if stat == JobStatus.terminated :
                    raise RuntimeError ("find terminated job. exit. should restart" )
                if stat != JobStatus.finished :
                    find_unfinish = True
                    break
            if find_unfinish == False :
                return
            else :
                time.sleep (60)
            print ("# checked at " + time.strftime("%Y-%m-%d %H:%M:%S"))            
    
    def submit_string (self,
                       step,
                       all_at_once = False,
                       ) :
        string_name = self.mk_string_name (step)
        if False == os.path.exists (string_name) :
            raise RuntimeError ("no string dir " + string_name + " should be generated before submitting")
        string = np.loadtxt (string_name + "/string.out")
        os.chdir (string_name)
        string_dir = os.getcwd() + "/"
        job_list = []
        job_id_list = []
        for node_idx in range (string.shape[0]) :            
            node_name = self.mk_node_name (node_idx)
            os.chdir (node_name)
            node_dir = os.getcwd() + "/"
            dep_idx = self.mk_dep_idx (node_idx)
            dep_finished = False
            if dep_idx >= 0 and job_list[dep_idx].check_status() == JobStatus.finished :                
                dep_finished = True
            print (node_name + " dep " + str(dep_idx) + " finished " + str(dep_finished))            
            if (all_at_once) or (step != 0 or dep_idx < 0) or dep_finished :
                sp.check_call ([string_dir+self.cmd_job_scpt])
            else :
                dep_job_id = job_id_list[dep_idx]
                sp.check_call ([string_dir+self.cmd_job_scpt, dep_job_id] )
            job = PBSJob (node_dir, "PBS.sub")
            job_list.append (job)
            job_id = job.submit ()
            job_id_list.append (job_id)
            os.chdir (string_dir)  
        os.chdir (self.base_dir)
        return job_list
        
    def generate_string (self,
                         step,
                         string
                         ) :
        os.chdir (self.base_dir)
        string_name = self.mk_string_name (step)
        if False == os.path.exists (string_name) :
            ret = Popen(["cp", '-a', self.string_template, string_name], stdout=PIPE, stderr=PIPE)
            stdout, stderr = ret.communicate()
            if ret.returncode != 0 :
                raise RuntimeError ("cannot copy template dir to " + string_name)
        np.savetxt (string_name + "/string.out", string)
        self.generate_string_node_dir (step, string)

        
    def generate_string_node_dir (self,
                                  step,
                                  string
                                  ) :
        string_name = self.mk_string_name (step)
        if False == os.path.exists (string_name) :
            raise RuntimeError ("the string path " + string_name + " should be made before making the nodes")
        os.chdir (string_name)
        if step == 0 : 
            for node_idx in range (string.shape[0]) :
                node_center = string[node_idx]
                self.mk_node_param (node_center)
                node_name = self.mk_node_name (node_idx)
                dep_idx = self.mk_dep_idx (node_idx)
                if dep_idx < 0 :
                    ret = Popen ([self.cmd_gen_dir, node_name],  stdout=PIPE, stderr=PIPE)
                else :
                    dep_node_name = self.mk_node_name (dep_idx)
                    ret = Popen ([self.cmd_gen_dir, node_name, dep_node_name],  stdout=PIPE, stderr=PIPE)
                stdout, stderr = ret.communicate()
                if ret.returncode == 1 :
                    raise RuntimeError ("cannot generate node. LOCATION: string: " +
                                        string_name +
                                        " node index: " +
                                        str(node_idx) +
                                        ".  error info: " +
                                        str(stderr, encoding='ascii') )
                if ret.returncode == 10 :
                    print ("# detected node: " + string_name + "/" + node_name + " updated.")
                    sp.check_call ([self.cmd_update_dir, node_name])
                else :
                    print ("# generated node: " + string_name + "/" + node_name)
        else :
            prev_string_name = self.mk_string_name (step-1)
            for node_idx in range (string.shape[0]) :
                node_center = string[node_idx]
                self.mk_node_param (node_center)
                node_name = self.mk_node_name (node_idx)
                dep_node_name = "../" + prev_string_name + "/" + node_name
                ret = Popen ([self.cmd_gen_dir, node_name, dep_node_name],  stdout=PIPE, stderr=PIPE)
                stdout, stderr = ret.communicate()
                if ret.returncode == 1 :
                    raise RuntimeError ("cannot generate node. LOCATION: string: " +
                                        string_name +
                                        " node index: " +
                                        str(node_idx) +
                                        ".  error info: " +
                                        str(stderr, encoding='ascii') )
                if ret.returncode == 10 :
                    print ("# detected node: " + string_name + "/" + node_name + " updated.")
                    sp.check_call ([self.cmd_update_dir, node_name])
                else :
                    print ("# generated node: " + string_name + "/" + node_name)
        os.chdir (self.base_dir)

    def mk_dep_idx (self,
                    idx,
                    ) :
        dep_idx = int(idx / self.dep_sec_size)
        dep_idx = dep_idx * self.dep_sec_size
        return dep_idx - 1

    def mk_node_param (self,
                       node_center
                       ) :        
        center_string = ""
        for idx in range (np.size(node_center)) :
            center = str(node_center[idx])
            if idx == 0 :
                center_string = center_string + str(center)
            else :
                center_string = center_string + ',' + str(center)                
        self.replace ("parameters.sh", "res_centers=.*", "res_centers=" + center_string)

    def mk_string_name (self,
                        step
                        ) :
        return "string.%06d" % (step)

    def mk_node_name (self,
                      step
                      ) :
        return "node.%06d" % (step)

    def replace (self, file_name, pattern, subst) :
        file_handel = open (file_name, 'r')
        file_string = file_handel.read ()
        file_handel.close ()
        file_string = ( re.sub (pattern, subst, file_string) )
        file_handel = open (file_name, 'w')
        file_handel.write (file_string)
        file_handel.close ()
        
if __name__ == "__main__" :
    node_start  = np.array([0.410,0.464,0.632])
    node_end    = np.array([0.012,0.028,0.021])
    numb_node   = 21
    dim = np.size(node_start)
    if dim != np.size(node_end):
        raise NameError ('The dimention of starting and ending nodes of the string should match!')
    string = np.zeros ((dim, numb_node))
    for ii in range (dim):
        string [ii] = np.linspace (node_start[ii], node_end[ii], numb_node)
    string = string.T
    print (str(string))
    sf = StringForce ("template.string")
    node = sf.compute (0, string)
    node = sf.compute (1, string)
    print (str(node))
