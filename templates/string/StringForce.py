#!/usr/bin/env python3

import os
import sys
import logging
Logger = logging.getLogger("string_method")
from subprocess         import Popen, PIPE
import re
import numpy as np
# import BatchJob

class StringForce (object) :
    def __init__ (self,
                  string_template
                  ) :
        self.base_dir = os.getcwd() + '/'
        self.string_template = string_template
        self.cmd_gen_dir = "tools/gen.dir.sh"
        self.cmd_job_scpt = "tools/mk.batch.sub.sh"

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
                 string) :
        self.set_up_string (step, string)
        self.do_simulation ()
        force = self.analyze ()
        return force

    def set_up_string_first (self,
                             string
                             ) :
        step = 0
        self.generate_string_dir (step)
        self.generate_string_node_dir (step, string)

    def submit_string (self,
                       step,
                       string
                       ) :
        string_name = self.mk_string_name (step)
        if False == os.path.exists (string_name) :
            Logger.error ("no string dir " + string_name + " should be generated before submitting")
            sys.exit(1)
            
        
        
    def generate_string_dir (self,
                             step
                             ) :
        string_name = self.mk_string_name (step)
        os.chdir (self.base_dir)
        if False == os.path.exists (string_name) :
            ret = Popen(["cp",'-a',self.string_template,string_name], stdout=PIPE, stderr=PIPE)
            stdout, stderr = ret.communicate()
            if ret.returncode != 0 :
                Logger.error ("cannot copy template dir to " + string_name)
                sys.exit (1)
        return string_name
    
    def generate_string_node_dir (self,
                                  step,
                                  string
                                  ) :
        string_name = self.mk_string_name (step)
        if False == os.path.exists (string_name) :
            generate_string_dir (step)
        os.chdir (string_name)
        if step == 0 :
            prev_node_dir = ""
            for node_idx in range (string.shape[0]) :
                node_center = string[node_idx]
                print ("index " + str(node_idx) + " center " + str(node_center))
                self.mk_node_param (node_idx, node_center)
                if node_idx == 0 :
                    ret = Popen ([self.cmd_gen_dir],  stdout=PIPE, stderr=PIPE)
                else :
                    ret = Popen ([self.cmd_gen_dir, prev_node_dir],  stdout=PIPE, stderr=PIPE)
                stdout, stderr = ret.communicate()
                if ret.returncode == 1 :
                    Logger.error ("cannot generate node. LOCATION: string: " + string_name + " node index: " + str(node_idx))
                    sys.exit(1)
                node_prefix = self.mk_node_name (node_idx)
                ret = Popen (["ls | grep " + node_prefix], stdout=PIPE, stderr=PIPE, shell = True)
                stdout, stderr = ret.communicate()
                if ret.returncode != 0 :
                    print (str(stdout, encoding='ascii'))                    
                    print (str(stderr, encoding='ascii'))
                    Logger.error ("cannot grep node name. LOCATION: string: " + string_name + " node index: " + str(node_idx))
                    sys.exit(1)
                prev_node_dir = str(stdout, encoding='ascii').split ('\n')[0]
                print ("# generated node: " + string_name + "/" + prev_node_dir)
        os.chdir (self.base_dir)

    def mk_node_param (self,
                       node_idx,
                       node_center
                       ) :        
        node_prefix = self.mk_node_name (node_idx)
        self.replace ("parameters.sh", "dir_prefix=.*", "dir_prefix=" + node_prefix)
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
    node_start  = np.array([0.012,0.028,0.021])
    node_end    = np.array([0.410,0.464,0.632])
    numb_node   = 3
    dim = np.size(node_start)
    if dim != np.size(node_end):
        raise NameError ('The dimention of starting and ending nodes of the string should match!')
    string = np.zeros ((dim, numb_node))
    for ii in range (dim):
        string [ii] = np.linspace (node_start[ii], node_end[ii], numb_node)
    string = string.T
    print (str(string))
    sf = StringForce ("restraint")
    node = sf.set_up_string_first (string)
    print (str(node))
