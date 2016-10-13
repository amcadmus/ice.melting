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

def select_node_from_list (gened_list,
                           ungened,
                           radius) :
    ### divide the ungened to select and unselect groups
    select = [[]]
    select_index = [[]]
    unselect = [[]]
    select_mark = np.zeros (ungened.shape[0]).astype (int)

    for ii in range (ungened.shape[0]) :
        for jj in range (len(gened_list)) :
            gened = gened_list[jj]
            assert (gened.shape[1] == ungened.shape[1])
            ret = dist (gened, ungened[ii])
            if ret[0] <= radius :
                if np.size(select) == 0 :
                    select = np.array([ungened[ii]])
                    select_index = np.array([[jj, ret[1]]])
                else :
                    select = np.append (select, [ungened[ii]], axis = 0)
                    select_index = np.append (select_index, np.array([[jj, ret[1]]]), axis = 0)
                select_mark[ii] = 1
        # print ("sel %d " % np.size(select))
        # print ("unsel %d" % unselect.shape[0])
    for ii in range (ungened.shape[0]) :
        if select_mark[ii] != 1 :
            if np.size(unselect) == 0 :
                unselect = [ungened[ii]]
            else :
                unselect = np.append (unselect, [ungened[ii]], axis = 0)
    return [select, select_index, unselect]

def make_dir_name (step) :
    return "step.%06d" % (step)

# def generate_dir (sel,
#                   dir_name,
#                   dep_dir_name) :
#     sf = StringForce ("template.string")
#     dep_nodes = np.loadtxt (dep_dir_name + "/string.out")

#     if not os.path.exists (dir_name) :
#         sp.check_call ("cp -a " + sf.string_template + " " + dir_name, shell = True)
#     np.savetxt (dir_name + "/string.out", sel)

#     base_dir = os.getcwd() + '/'
#     os.chdir (dir_name)
#     for ii in range (sel.shape[0]) :
#         node_center = sel[ii]
#         ret = dist (dep_nodes, node_center)
#         dep_posi = ret[1]
#         dep_node_name = sf.mk_node_name (dep_posi)
#         dep_node_name = "../" + dep_dir_name + "/" + dep_node_name
#         node_name = sf.mk_node_name (ii)
#         sf.mk_node_param (node_center)
#         ret = Popen ([sf.cmd_gen_dir, node_name, dep_node_name],  stdout=PIPE, stderr=PIPE)
#         stdout, stderr = ret.communicate()
#         if ret.returncode == 1 :
#             raise RuntimeError ("cannot generate node. LOCATION: string: " +
#                                 dir_name +
#                                 " node index: " +
#                                 str(ii) +
#                                 ".  error info: " +
#                                 str(stderr, encoding='ascii') +
#                                 ".  error info: " +
#                                 str(stdout, encoding='ascii')
#                                 )
#         if ret.returncode == 10 :
#             logging.info ("# detected  node: " + dir_name + "/" + node_name + " do nothing.")
#         else :
#             logging.info ("# generated node: " + dir_name + "/" + node_name + " with dep " + dep_node_name)
            
#     os.chdir (base_dir)
                  
def generate_dir (select,
                  select_index,
                  gened_list) :
    sf = StringForce ("template.string")

    step = len(gened_list)
    dir_name = make_dir_name (step)
    if not os.path.exists (dir_name) :
        sp.check_call ("cp -a " + sf.string_template + " " + dir_name, shell = True)
    np.savetxt (dir_name + "/string.out", select)
    
    base_dir = os.getcwd() + '/'
    os.chdir (dir_name)
    for ii in range (select.shape[0]) :
        dep_step_index = select_index[ii][0]
        dep_node_index = select_index[ii][1]
        pair_dist = np.linalg.norm (select[ii] - np.array(gened_list[dep_step_index][dep_node_index]))
        dep_node_name = "../" + make_dir_name(dep_step_index) + "/" + sf.mk_node_name(dep_node_index)
        node_name = sf.mk_node_name (ii)
        sf.mk_node_param (select[ii])
        ret = Popen ([sf.cmd_gen_dir, node_name, dep_node_name],  stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if ret.returncode == 1 :
            raise RuntimeError ("cannot generate node. LOCATION: string: " +
                                dir_name +
                                " node index: " +
                                str(ii) +
                                ".  error info: " +
                                str(stderr, encoding='ascii') +
                                ".  error info: " +
                                str(stdout, encoding='ascii')
                                )
        if ret.returncode == 10 :
            logging.info ("# detected  node: " + dir_name + "/" + node_name + " do nothing.")
        else :
            logging.info ("# generated node: " + dir_name + "/" + node_name + " with dep " + dep_node_name
                          + " dist " + str(pair_dist))
            
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
        
        select = np.loadtxt (make_dir_name(0)  + "/string.out")
        gened = [select]
        unselect = np.loadtxt (args.input)
        radius = max_spacing * 1.1
        step = 1
        while True :
            ret = select_node_from_list (gened, unselect, radius)
            select = ret[0]
            select_index = ret[1]
            unselect = ret[2]            
            dir_name = make_dir_name (step)            
            print ("step %d" % step)
            print (select.shape)
            print (select_index.shape)
            generate_dir (select, select_index, gened)
            gened.append (select)            
            logging.info ("numb generated list %d", len(gened))
            if np.size(unselect) == 0 :
                break
            step = step + 1
            
        

if __name__ == "__main__":

    main ()
