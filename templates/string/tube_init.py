#!/usr/bin/env python3

import math
import os
import sys
import argparse
import numpy as np
import scipy 
import logging
from scipy.interpolate import interp1d
from scipy.optimize import minimize

def unique_rows (a) :
    if a.shape[1] <= 1 :
        return a
    a = np.ascontiguousarray(a)
    unique_a = np.unique (a.view([('', a.dtype)] * a.shape[1]))
    return unique_a.view (a.dtype).reshape(unique_a.shape[0], a.shape[1])

def interp_string (string) :    
    numb_node = string.shape[0]
    alpha     = np.linspace (0, 1, numb_node)
    alpha_seg = np.linspace (0, 1, numb_node)
    alpha_seg[0] = 0
    for jj in range (1, numb_node):
        alpha_seg[jj] = np.linalg.norm (string[jj] - string[jj-1])
    alpha = np.cumsum (alpha_seg)
    alpha = alpha / alpha[-1]
    return interp1d (alpha, string, axis=0, kind="cubic")

class StringPointPointDist (object) :
    def __init__ (self, 
                  string,
                  point ) :
        # self.string = string
        self.point = point
        # self.smooth_string = interp_string (string)
        self.smooth_string = string

    def dist (self, 
              alpha) :
        str_point = self.smooth_string (alpha)
        return np.linalg.norm (str_point - self.point)
        

class StringPointDist (object) :
    def __init__ (self,
                  string, 
                  fixed_point,
                  fixed_dim = []) :
        self.string = string
        self.fixed_point = fixed_point
        self.fixed_dim = fixed_dim        
        self.numb_dim = len(string(0))
        self.numb_var_dim = self.numb_dim - len(self.fixed_dim)
        self.var_dim = []
        for ii in range(self.numb_dim) :
            if not ii in self.fixed_dim :
                self.var_dim.append (ii)

    def dist_full_point (self, 
                         point) :
        spd = StringPointPointDist (self.string, point)
        res = scipy.optimize.minimize_scalar (spd.dist, bounds = (0, 1), method = 'bounded')
        # print (res)
        # for ii in range(100) :
        #     hh = 1./100.
        #     print (np.linalg.norm(spd.smooth_string(hh * ii) - point))
        # print (string[0])
        # print (point)
        # print ("%f %f" % (spd.dist(0), np.linalg.norm(string[0] - point)) ) 
        # print ("%f %f" % (spd.dist(1), np.linalg.norm(string[-1] - point)) )
        # print ("%f %f" % (spd.dist(0.4444444444444), np.linalg.norm(string[4] - point)) )
        return [res.fun, res.x]

    def dist (self,
              point) :
        assert (len(point) == self.numb_var_dim)
        string_point = np.array (self.fixed_point, copy = True) 
        for ii in range (self.numb_var_dim) :
            string_point[self.var_dim[ii]] = point[ii]
        # print (self.fixed_point)
        # print (self.dist_full_point (self.fixed_point))
        # print (string_point)
        # print (self.dist_full_point (string_point))
        return self.dist_full_point (string_point)[0]

def gen_point_var_dim (string, 
                       point,
                       fixed_dim = []) :
    myspd = StringPointDist (string, point, fixed_dim)
    start_point = []
    for ii in range (myspd.numb_var_dim) :
        start_point.append (point[myspd.var_dim[ii]])
    # print (point)
    # print (start_point)
    # print (myspd.dist (start_point))    
    res = scipy.optimize.minimize (myspd.dist, start_point, method = 'nelder-mead',options={'xtol': 1e-8})
    # return [res.fun, res.x]
    # print (res.fun)
    retpnt = np.array (point, copy = True)
    for ii in range (myspd.numb_var_dim) : 
        retpnt[myspd.var_dim[ii]] = res.x[ii]
    return retpnt

def gen_mesh_node_index (string_node,
                         mesh_spacing,
                         tube_radius, 
                         constrained_dim = []) :
    numb_dim = string_node.shape[0]
    assert (numb_dim == mesh_spacing.shape[0])
    mesh_invh = 1./mesh_spacing
    index_string_node = string_node * mesh_invh
    index_string_node = index_string_node.astype (int)
    print ("gen mesh for string node %s" % string_node)
    range_lower = np.zeros (string_node.shape).astype(int)
    range_upper = np.ones  (string_node.shape).astype(int)
    for ii in range (numb_dim) :
        if ii in constrained_dim : 
            continue
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
            if dd in constrained_dim :
                node_index[dd] = -1
            else :
                node_index[dd] = value % range_size[dd] + range_lower[dd]
            value = value / range_size[dd]
#        print (node_index)
        if np.size(mesh_node) == 0 :
            mesh_node = np.array([node_index])
        else :
            mesh_node = np.concatenate ((mesh_node, [node_index]), axis = 0)

    # print (mesh_node)
    # print (mesh_node.shape)
    return mesh_node

def gen_tube (string,
              mesh_spacing,
              tube_radius, 
              constrained_dim = []) :
    numb_node = string.shape[0]
    numb_dims = string.shape[1]
    tube_index = np.array([[]])
#    for ii in range (2) :
    for ii in range (numb_node) :
        mesh_node = gen_mesh_node_index (string[ii], mesh_spacing, tube_radius, constrained_dim)
#        print (mesh_node)
        if np.size(tube_index) == 0 :
            tube_index = mesh_node
        else :
            tube_index = np.concatenate ((tube_index, mesh_node), axis = 0)
        tube_index = unique_rows (tube_index)

    fixed_dim = []
    for ii in range(numb_dims) :
        if not ii in constrained_dim :
            fixed_dim.append (ii)
    # print (fixed_dim)

    smooth_string = interp_string (string)

    tube = tube_index.astype (float)
    for ii in range(tube.shape[0]) :
        tmp_tube = tube[ii] * mesh_spacing
        if len(fixed_dim) == numb_dims :
            tube[ii] = tmp_tube
        else :
            for jj in constrained_dim :
                tmp_tube[jj] = smooth_string(0.5)[jj]
            tube[ii] = gen_point_var_dim (smooth_string, tmp_tube, fixed_dim)
            print ("%s %s" % (tmp_tube, tube[ii]))
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
    parser.add_argument('-c', '--constrained-dim', type=int, nargs = '*', default = [],
                        help='The constrained dimension')
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
    tube = gen_tube (string, mesh_spacing, args.tube_radius, args.constrained_dim)
    logging.info ("tube generated. numb_node %d " % tube.shape[0])
    print ("tube generated. numb_node %d " % tube.shape[0])

    # save tube
    fp = open (args.output, "w")
    for ii in args.constrained_dim :
        mesh_spacing[ii] = -1.
    fp.write ("# mesh_spacing %s\n" % mesh_spacing)
    fp.close()
    fp = open (args.output, "ab")
    np.savetxt (fp, tube)
    fp.close ()
    logging.info ("tube saved")

    # point = np.array([0.1, 0.2, 0.5, 0.5])
    # myspd = StringPointDist (string, point, fixed_dim = [1,3] )
    # mydist = myspd.dist_full_point (point)
    # print (mydist)
    # print ("test string point class")
    # point = np.array([0.3, 0.1])
    # myspd.dist (point)

    # maxval = 100.
    # x0 = np.linspace (0.1, 0.5, 100)
    # x2 = np.linspace (0.1, 0.5, 100)
    # for ii in x0 :
    #     for jj in x2 :
    #         mypt = [ii, point[1], jj, point[3]]
    #         # myval = myspd.dist (np.array([ii, jj]))
    #         myval = myspd.dist_full_point (mypt)[0]
    #         if (myval < maxval) :
    #             print ("%.3f %.3f  %f" % (ii, jj, myval) )
    #             maxval = myval
    # res = gen_point_var_dim (string, point, fixed_dim = [1,3])
    # print (point)
    # print (res)
    # print (myspd.dist_full_point(res)[0])

if __name__ == '__main__':
    main()
    
