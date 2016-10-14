#!/usr/bin/env python3

import math
import os
import sys
import argparse
import numpy as np
import logging

def unique_rows (a) :
    if a.shape[1] <= 1 :
        return a
    a = np.ascontiguousarray(a)
    unique_a = np.unique (a.view([('', a.dtype)] * a.shape[1]))
    return unique_a.view (a.dtype).reshape(unique_a.shape[0], a.shape[1])

def gen_mesh_node_index (string_node,
                         mesh_spacing,
                         tube_radius) :
    numb_dim = string_node.shape[0]
    assert (numb_dim == mesh_spacing.shape[0])
    mesh_invh = 1./mesh_spacing
    index_string_node = string_node * mesh_invh
    index_string_node = index_string_node.astype (int)
    print ("gen mesh for string node %s" % string_node)
    range_lower = np.zeros (string_node.shape).astype(int)
    range_upper = np.zeros (string_node.shape).astype(int)
    for ii in range (numb_dim) :
        index_radius = int(tube_radius / mesh_spacing[ii])
        if index_radius < 1 :
            index_radius = 1
        range_lower[ii] = index_string_node[ii] - index_radius
        range_upper[ii] = index_string_node[ii] + index_radius + 2
        if range_lower[ii] < 0 :
            range_lower[ii] = 0
        if range_upper[ii] < 0 :
            range_upper[ii] = 0
    range_size = range_upper - range_lower
    numb_node = np.cumprod(range_size)[-1]    
#    print (range_size)
    mesh_node = [[]]
    for ii in range (numb_node) :
        node_index = np.zeros (numb_dim).astype (int)
        value = ii
        for dd in range (numb_dim) :
            node_index[dd] = value % range_size[dd] + range_lower[dd]
            value = value / range_size[dd]
#        print (node_index)
        if np.size(mesh_node) == 0 :
            mesh_node = np.array([node_index])
        else :
            mesh_node = np.concatenate ((mesh_node, [node_index]), axis = 0)

#    print (mesh_node)
    return mesh_node

def gen_tube (string,
              mesh_spacing,
              tube_radius) :
    numb_node = string.shape[0]
    numb_dims = string.shape[1]
    tube_index = np.array([[]])
#    for ii in range (2) :
    for ii in range (numb_node) :
        mesh_node = gen_mesh_node_index (string[ii], mesh_spacing, tube_radius)
#        print (mesh_node)
        if np.size(tube_index) == 0 :
            tube_index = mesh_node
        else :
            tube_index = np.concatenate ((tube_index, mesh_node), axis = 0)
        tube_index = unique_rows (tube_index)

    tube = tube_index.astype (float)
    for ii in range(tube.shape[0]) :
        tube[ii] = tube[ii] * mesh_spacing
    return tube

def main () :
    parser = argparse.ArgumentParser(
        description="*** Initialize a tube around a string. Do not do the simulation. ***")

    parser.add_argument('-s', '--string', default = "step.000000",
                        help='The string. ')
    parser.add_argument('-o', '--output', default = "tube.out",
                        help='The file for output tube. ')
    parser.add_argument('-d', '--mesh-spacing', type=float, nargs = '*', default = 0.002,
                        help='The mesh spacing')
    parser.add_argument('-r', '--tube-radius', type=float, default = 0.02,
                        help='The radius of the tube around the string')
    args = parser.parse_args()

    logging.basicConfig (filename="tube_init.log", filemode="w", level=logging.INFO, format='%(asctime)s %(message)s')

    # parse the string
    logging.info ("gen tube for string: %s" % args.string)
    string = np.loadtxt (args.string + "/string.out")
    numb_node = string.shape[0]
    numb_dims = string.shape[1]
    logging.info ("numb_node of the string: %s" % numb_node)
    logging.info ("numb_dims of the string: %s" % numb_dims)

    # mesh spacing
    if np.size(args.mesh_spacing) != 1 and numb_dims != np.size(args.mesh_spacing) :
        logging.fatal ("unmatched string dim and the dim in mesh size, die")
        raise RuntimeError ("unmatched string dim and the dim in mesh size, die")
    if np.size(args.mesh_spacing) == 1 :
        mesh_spacing = np.repeat (args.mesh_spacing, numb_dims)
    else :
        mesh_spacing = np.array (args.mesh_spacing)
    logging.info ("mesh spacing %s" % mesh_spacing)
    logging.info ("tube radius  %s" % args.tube_radius)

    # gen tube
    tube = gen_tube (string, mesh_spacing, args.tube_radius)
    logging.info ("tube generated. numb_node %d " % tube.shape[0])
    print ("tube generated. numb_node %d " % tube.shape[0])

    # save tube
    fp = open (args.output, "w")
    fp.write ("# mesh_spacing %s\n" % mesh_spacing)
    fp.close()
    fp = open (args.output, "ab")
    np.savetxt (fp, tube)
    fp.close ()
    logging.info ("tube saved")

if __name__ == '__main__':
    main()
