#!/usr/bin/env python3

import os
import argparse
import re
import numpy as np
import logging
import subprocess as sp
import enum as Enum
from StringForce import StringForce

def dist (gened,
          node ) :
    min_p = 0
    min_v = np.linalg.norm (gened[min_p] - node)
    for ii in range (1, gened.shape[0]) :
        new_v = np.linalg.norm (gened[ii] - node)
        if new_v < min_v :
            min_v = new_v
            min_p = ii
    return [min_v, min_p]

def select_node (gened,
                 ungened,
                 radius) :
    ### divide the ungened to select and unselect groups
    assert (gened.shape[1] == ungened.shape[1])
    select = [[]]
    unselect = [[]]
    for ii in range (ungened.shape[0]) :
        ret = dist (gened, ungened[ii])        
        if ret[0] <= radius :
            if np.size(select) == 0 :
                select = np.array([ungened[ii]])
            else :
                select = np.append (select, [ungened[ii]], axis = 0)
        else :
            if np.size(unselect) == 0 :
                unselect = np.array([ungened[ii]])
            else :
                unselect = np.append (unselect, [ungened[ii]], axis = 0)
        # print ("sel %d " % np.size(select))
        # print ("unsel %d" % unselect.shape[0])
    return [select, unselect]

def make_dir_name (step) :
    return "tubestep.%06d" % (step)

def generate_dir (sel,
                  dir_name,
                  dep_dir_name) :
    sf = StringForce ("template.string")
    dep_nodes = np.loadtxt (dep_dir_name + "nodes.out")

    base_dir = os.getcwd() + '/'
    os.chdir (dir_name)
    np.savetxt ("nodes.out", sel)

    for ii in range (sel.shape[0]) :
        node_center = sel[ii]
        ret = dist (dep_nodes, node_center)
        dep_posi = ret[1]
        dep_node_name = sf.mk_node_name (dep_posi)
        dep_node_name = "../" + dep_dir_name + "/" + dep_node_name
        node_name = df.mk_node_name (ii)
        sf.mk_node_param (node_center)
        ret = Popen ([sf.cmd_gen_dir, node_name, dep_node_name],  stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if ret.returncode == 1 :
            raise RuntimeError ("cannot generate node. LOCATION: string: " +
                                dir_name +
                                " node index: " +
                                str(ii) +
                                ".  error info: " +
                                str(stderr, encoding='ascii') )
        if ret.returncode == 10 :
            logging.info ("# detected  node: " + dir_name + "/" + node_name + " do nothing.")
        else :
            logging.info ("# generated node: " + dir_name + "/" + node_name)
            
    os.chdir (base_dir)
                  

def main ():
    parser = argparse.ArgumentParser(
        description="*** Sample a tube. ***")
    parser.add_argument('-s', '--string', default = "string.000000",
                        help='The string. ')
    parser.add_argument('-i', '--input', default = "tube.out",
                        help='The tube. ')
    parser.add_argument('-j', '--gen-jobs', action = "store_true", default = False,
                        help='Generate jobs for a tube. ')
    # parser.add_argument('--gen-dir-script', default = "tools/gen.dir.sh",
    #                     help='Generate script for job. ')
    parser.add_argument('-t','--md-time', type=int, default=20,
                        help='Physical time of MD simulation in unit of ps.')
    args = parser.parse_args()

    logging.basicConfig (filename="sample_tube.log", filemode="w", level=logging.INFO, format='%(asctime)s %(message)s')

    sf = StringForce ("template.string")
    sf.replace ("template.string/parameters.sh", "md_time=.*", "md_time=" + str(args.md_time))

    tag_gen_jobs = 'tag_gen_jobs'
    if args.gen_jobs :
        if os.path.isfile (tag_gen_jobs) :
            msg = "jobs already generated, if want to regenerate, remove the tag: " + tag_gen_jobs
            logging.error (msg)
            raise RuntimeError (msg)
        if not os.path.isfile (args.string + "/string.out") :
            raise RuntimeError ("no file " + args.string)        
        if not os.path.isfile (args.input) :
            raise RuntimeError ("no file " + args.input)        
        fp = open (args.input)
        head = fp.readline ()
        head = re.sub (r"\#|[a-z]|_|\[|\]|\n","",head)
        mesh_spacing = np.fromstring (head, dtype=float, sep=" ")
        print ("mesh spacing %s" % mesh_spacing)
        max_spacing = np.amax (mesh_spacing)
        print ("max mesh spacing %s" % max_spacing)

        gened = np.loadtxt (args.string + "/string.out")
        ungened = np.loadtxt (args.input)
        print (gened.shape )
        print (ungened.shape )
        radius = max_spacing * 1.1
        step = 0
        while True :
            ret = select_node (gened, ungened, radius)
            sel = ret[0]
            unsel = ret[1]
            print (step)
            if sel.shape[0] == 0 :
                break
            dir_name = make_dir_name (step)
            if step == 0:
                dep_dir_name = string
            else :
                dep_dir_name = make_dir_name (step-1)
            generate_dir (sel, dir_name, dep_dir_name)
            break
            step = step + 1
            
        

if __name__ == "__main__":
    sys.path.append ('/usr/lib/python2.7/site-packages')
    main ()
