#!/usr/bin/env python3

import math
import os
import sys
import argparse
import numpy as np
import string_utils
from StringForce import StringForce
from subprocess import Popen, PIPE
from scipy.interpolate import interp1d

def __available_method () :
    LIST = ["linear", "resample", "resample-reinit"]
    names=""
    for ii in LIST :
        names += ii + ", "
    return names

def init_linear_string (node_start, node_end, numb_node):
    """ init a linear string between the start and end nodes. """
    dim = np.size(node_start)
    if dim != np.size(node_end):
        raise NameError ('The dimention of starting and ending nodes of the string should match!')
    string = np.zeros ((dim, numb_node))
    for ii in range (dim):
        string [ii] = np.linspace (node_start[ii], node_end[ii], numb_node)
    return string.T

def init_source_string (string_dir, numb_node_tgt) :
    """ init a string from existing nodes"""
    if not os.path.isdir (string_dir) :
        raise RuntimeError ("Dir " + string_dir + " not found")
    file_name = string_dir + "/string.out"
    if not os.path.exists (file_name) :
        raise RuntimeError ("cannot find file " + file_name)

    string = np.loadtxt (file_name)
    alpha = string_utils.arc_norm (string)
    smooth_str = interp1d (alpha, string, axis=0, kind="linear")

    alpha_eq = np.linspace (0, 1, numb_node_tgt)
    string = smooth_str (alpha_eq)
    
    return string

def generate_from_source (source_string_dir_,
                          string) :
    # check the dir and string file
    if not os.path.isdir (source_string_dir_) :
        raise RuntimeError ("cannot find dir " + source_string_dir_)
    cwd = os.getcwd()
    os.chdir (source_string_dir_)
    source_string_dir = os.getcwd ()
    os.chdir (cwd)
    file_name = source_string_dir + "/string.out"
    if not os.path.exists (file_name) :
        raise RuntimeError ("cannot find file " + file_name)
    # load source string
    source_string = np.loadtxt (file_name)
    if source_string.shape[1] != string.shape[1] :
        raise RuntimeError ("size of the string nodes does not match")
    # mk dir of the string
    str_force = StringForce ()
    string_name = str_force.mk_string_name (0)
    if not os.path.isdir(string_name) :        
        ret = Popen(["cp", '-a', str_force.string_template, string_name], stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if ret.returncode != 0 :
            raise RuntimeError ("cannot copy template dir to " + string_name)
    # generate nodes
    cwd = os.getcwd()
    cmd_gen_dir = "tools/gen.dir.absdep.sh"
    os.chdir (string_name)

    my_alpha = np.linspace (0, 1, string.shape[0])
    source_alpha = np.linspace (0, 1, source_string.shape[0])
    for ii in range (string.shape[0]) :
        # min_val = 1e10
        # min_posi = 0
        # for jj in range(source_string.shape[0]) :
        #     norm = np.linalg.norm (string[ii] - source_string[jj])
        #     # print ("id " + str(jj) +
        #     #        " diff " + str(string[ii]) + " " + str(source_string[jj]) + " norm " + str(norm))
        #     if norm < min_val :
        #         min_val = norm
        #         min_posi = jj
        for jj in range(source_string.shape[0]) :
            if (my_alpha[ii] >= source_alpha[jj]) :
                min_posi = jj
        this_node = str_force.mk_node_name (ii)
        str_force.mk_node_param (string[ii])
        source_node = str_force.mk_node_name (min_posi)
        ret = Popen ([cmd_gen_dir, this_node, source_string_dir+"/"+source_node],  stdout=PIPE, stderr=PIPE)
        stdout, stderr = ret.communicate()
        if ret.returncode == 1 :
            raise RuntimeError ("cannot generate node. LOCATION: string: " +
                                string_name +
                                " node index: " +
                                str(this_node) +
                                ".  error info: " +
                                str(stderr, encoding='ascii') )
        if ret.returncode == 10 :
            print ("# detected node: " + string_name + "/" + this_node)
        else :
            print ("# generated node: " + string_name + "/" + this_node)
    np.savetxt ("string.out", string)
    os.chdir (cwd)

def main ():
    parser = argparse.ArgumentParser(
        description="*** Initialize a string. ***")

    parser.add_argument('-m', '--method', default = "",
                        help='Method for generating the string. Available methods are: ' + __available_method())
    parser.add_argument('-n', '--numb-nodes', type=int, default = 20,
                        help='Number of nodes on a string.')
    parser.add_argument('-t','--md-time', type=int, default=20,
                        help='Physical time of MD simulation in unit of ps.')
    parser.add_argument('-d','--dep-size', type=int, default=1,
                        help='Number of nodes that depends on one node.')

    lg = parser.add_argument_group ("Linear string")
    lg.add_argument('-b', '--begin', type=float, nargs = '*',
                    help='Start of the string.')
    lg.add_argument('-e', '--end', type=float, nargs = '*',
                    help='End of the string.')

    sg = parser.add_argument_group ("Source string")
    sg.add_argument('-s', '--source', 
                    help='Generate a new string from this string.')

    args = parser.parse_args()

    str_force = StringForce ("template.string", args.dep_size)
    str_force.replace ("template.string/parameters.sh", "md_time=.*", "md_time=" + str(args.md_time))

    if args.method == "linear" :
        if args.begin == None :
            raise RuntimeError ("Begin of the string is empty")
        if args.end == None :
            raise RuntimeError ("End of the string is empty")
        string = init_linear_string (args.begin, args.end, args.numb_nodes)
        print (str(string))    
        str_force.generate_string (0, string)
        job = str_force.submit_string (0, False)
    elif args.method == "resample" :
        if args.source == None :
            raise RuntimeError ("No source string")
        string = init_source_string (args.source, args.numb_nodes)
        generate_from_source (args.source, string)
        job = str_force.submit_string (0, True)
    elif args.method == "resample-reinit" :
        if args.source == None :
            raise RuntimeError ("No source string")
        string = init_source_string (args.source, args.numb_nodes)
        str_force.generate_string (0, string)
        job = str_force.submit_string (0, False)        
    else :
        raise RuntimeError ("unknow method to generate the string!")

    str_force.wait_string (job)
    str_force.write_tag (0)
    force = str_force.statistic_string (0)
        
if __name__ == '__main__':
    main()
