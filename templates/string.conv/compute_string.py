#!/usr/bin/env python3
# to be converted to python3

import os
import argparse
import numpy as np
import subprocess as sp
import logging
from scipy.interpolate import interp1d
from subprocess import Popen, PIPE


def test_compute_force (xx):
    yy = np.zeros (np.size(xx))
    yy[0] = 8*xx[0]*(1-xx[0]*xx[0]-xx[1]*0.5)
    yy[1] = 2*(1-xx[0]*xx[0]-xx[1])
    return yy

class TestRandomForce (object) :
    """
    The test_compute_force + random noise
    """
    def __init__ (self,
                  sigma) :
        self.sigma = sigma

    def compute_force (self,
                       step, 
                       string ) :
        numb_node = string.shape[0]
        dim = string.shape[1]
        force = np.zeros ((numb_node, dim))
        for jj in range (numb_node):
            force[jj] = test_compute_force (string[jj])
        noise = np.random.normal (0, self.sigma, (numb_node, dim))
        np.savetxt ("string/force.%06d.out" % step, force + noise)
        return force + noise

def init_linear_string (node_start, node_end, numb_node):
    """ init a linear string between the start and end nodes. """
    dim = np.size(node_start)
    if dim != np.size(node_end):
        raise NameError ('The dimention of starting and ending nodes of the string should match!')
    string = np.zeros ((dim, numb_node))
    for ii in range (dim):
        string [ii] = np.linspace (node_start[ii], node_end[ii], numb_node)
    return string.T

def compute_string_tegent (alpha,
                           string,
                           delta_a = 0.001
                           ):
    """ compute the tangent vector of the string, it is normalized """
    tangent = np.zeros (string.shape)
    numb_node = string.shape[0]
    dim = string.shape[1]
    smooth_str = interp1d (alpha, string, axis=0, kind="cubic")
    tangent[0]  = ( smooth_str(alpha[ 0] + delta_a) - smooth_str(alpha[ 0]) ) / delta_a
    tangent[-1] = ( smooth_str(alpha[-1]) - smooth_str(alpha[-1] - delta_a) ) / delta_a
    for ii in range (1, numb_node-1):
        tangent[ii] = ( smooth_str(alpha[ii] + delta_a) - smooth_str(alpha[ii] - delta_a ) ) / (delta_a * 2.)
    norm_t = np.sqrt (np.sum(np.multiply (tangent, tangent), axis=1))
    for ii in range (numb_node):
        tangent[ii] = tangent[ii] / norm_t[ii]
    return tangent

def string_update_rhs (compute_force,
                       dt,
                       step,
                       string):
    """ compute the dt * force """
    return dt * compute_force (step, string)

def update_string_Euler (compute_force,
                         dt,
                         step,
                         string):
    incr = string_update_rhs (compute_force, dt, step, string)
    return string + incr

def update_string_RK2 (compute_force, dt, step, string):
    my_step = int(step * 2)
    in_k1 = string
    k1 = string_update_rhs (compute_force, dt, my_step+0, in_k1)
    in_k2 = string + 0.5 * k1
    k2 = string_update_rhs (compute_force, dt, my_step+1, in_k2)
    return string + k2

def update_string_RK4 (compute_force, dt, step, string):
    my_step = int(step * 4)
    in_k1 = string
    k1 = string_update_rhs (compute_force, dt, my_step+0, in_k1)
    in_k2 = string + 0.5 * k1
    k2 = string_update_rhs (compute_force, dt, my_step+1, in_k2)
    in_k3 = string + 0.5 * k2
    k3 = string_update_rhs (compute_force, dt, my_step+2, in_k3)
    in_k4 = string + 1.0 * k3
    k4 = string_update_rhs (compute_force, dt, my_step+3, in_k4)
    return string + (1./6.) * k1 + (1./3.) * k2 + (1./3.) * k3 + (1./6.) * k4

def compute_string (compute_force,              # function for computing the force
                    string,                     # the input string
                    dt = 0.05,                  # artificial time step for updating the string
                    max_iter = 200,             # maximum allowed number of iterations
                    start_iter = 0
                    ):            
    """ compute the string"""
    sp.check_call ("rm -fr string/string*out", shell = True)
    factor_Q = 1.5
    numb_node = string.shape[0]
    dim = np.size(string[0])
    # check validity of the inputs
    if dim != np.size(string[-1]):
        raise NameError ('The dimention of starting and ending nodes of the string should match!')
    if numb_node <= 2:
        raise NameError ('The number of nodes on string should be larger than 2')
    # initialize
    alpha       = np.linspace (0, 1, numb_node)
    alpha_seg   = np.linspace (0, 1, numb_node)
    alpha_eq    = np.linspace (0, 1, numb_node)
    incr_hist   = [[]]

    conv_file = open ("conv.out", "w")
    # starts the main loop
    for ii in range (start_iter, max_iter):
        # update the string
        string = update_string_Euler (compute_force, dt, ii, string)
        # string = update_string_RK4 (compute_force, dt, string)
        # update the arc
        alpha_seg[0] = 0
        for jj in range (1, numb_node):
            alpha_seg[jj] = np.linalg.norm (string[jj] - string[jj-1])
        alpha = np.cumsum (alpha_seg)
        alpha = alpha / alpha[-1]
        # reparameterize the string if needed
        smooth_str = interp1d (alpha, string, axis=0, kind="cubic")
        if np.max (alpha_seg[1:]) / np.min(alpha_seg[1:]) > factor_Q :
            string = smooth_str (alpha_eq)
            logging.info ("string %06d: resampled .", ii+1)
        # print string
        np.savetxt ("string/string.%06d.out" % ii, string)
        # compute the max norm force as measure of convergence
        if ii != start_iter :
            norm_string = smooth_str (alpha_eq)
            diff_string = norm_string - norm_string_old
            norm_string_old = np.copy (norm_string)
            diff = np.sqrt (np.sum (np.multiply (diff_string, diff_string), axis=1))                        
            diff_inf = np.max( diff )
            new_item = np.array([ii, diff_inf])
            new_item = new_item[np.newaxis,:]
            if np.size (incr_hist) == 0:
                incr_hist = new_item
            else:
                incr_hist = np.append (incr_hist, new_item, axis=0)
            logging.info ("string %06d: updated with timestep %e . String difference is %e", ii+1, dt, diff_inf)
            conv_file.write (str(ii) + " " + str(diff_inf) + "\n")
        else :
            norm_string_old = smooth_str (alpha_eq)            
            logging.info ("string %06d: updated with timestep %e .", ii+1, dt)
#    print incr_hist
    conv_file.close ()
    return string    

def main ():

    parser = argparse.ArgumentParser(
        description="*** Initialize a string. ***")

    parser.add_argument('-n', '--nodes', type=int, default = 200,
                        help='The number of nodes on the string')
    parser.add_argument('-s', '--sigma', type=float, default = 1e-3,
                        help='The standard deviation of the noise on force')
    parser.add_argument('-S', '--seed', type=int, default = 0,
                        help='The seed for random generator')
    parser.add_argument('-d', '--step', type=float, default = 2e-6,
                        help='Step size for evolving string')
    parser.add_argument('-m', '--max-step', type=int, default = 300,
                        help='Maximum number of steps')
    parser.add_argument('-t','--md-time', type=int, default=20,
                        help='Physical time of MD simulation in unit of ps.')

    args = parser.parse_args()

    logging.basicConfig (filename="compute_string.log", filemode="a", level=logging.INFO, format='%(asctime)s %(message)s')

    np.random.seed (args.seed)
    string = init_linear_string ([-1,0], [1,0], args.nodes)
    force_noise = TestRandomForce (args.sigma)
    string = compute_string (force_noise.compute_force, string, args.step, args.max_step, 0)
    
if __name__ == "__main__":
    main ()
