#!/usr/bin/env python3

import os
import logging
import sys
from subprocess         import Popen, PIPE
import BatchJob

class StringForce (object) :
    def __init__ (self,
                  string_template) :
        self.string_template = string_template
        if False == os.path.exists (self.string_template) :
            Logger.error ("cannot find string template " + self.string_template)
        os.chdir (self.string_template)
        if False == os.path.exists (parameters.sh) :
            Logger.error ("no file parameters.sh")
            sys.exit ()
        if False == os.path.exists (tools/gen.dir.sh) :
            Logger.error ("no file tools/gen.dir.sh")
            sys.exit ()
        if False == os.path.exists (tools/mk.batch.sub.sh) :
            Logger.error ("no file tools/mk.batch.sub.sh")
            sys.exit ()
    
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
        node_name = 
                 
                 
