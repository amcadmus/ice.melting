#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import bisect

global_line = 0
global_mapping = None

def arc_seg (string) :
    ### compute arc (alpha) segments
    numb_node = string.shape[0]
    alpha_seg = np.zeros (numb_node)
    alpha_seg[0] = 0
    for jj in range (1, numb_node):
        alpha_seg[jj] = np.linalg.norm (string[jj] - string[jj-1])
    return alpha_seg

def arc (string) :
    ### compute arc parameter
    alpha_seg = arc_seg (string)
    alpha = np.cumsum (alpha_seg)
    return alpha

def arc_norm (string) :
    ### compute normalized arc parameter
    alpha = arc (string)
    alpha = alpha / alpha[-1]
    return alpha

def __solve_func (x) :
    return global_mapping(x) - global_line

def resample_string (string_,
                     new_numb_node,
                     weighting_ = [[0, 1], [1, 1]]) :
    string = string_
    global global_line
    global global_mapping
    
    # build the smth string
    alpha = arc_norm (string)
    smooth_str = interp1d (alpha, string, axis=0, kind="linear")

    # smooth weight
    weighting = np.array(weighting_)
    smt_w = interp1d (weighting[:,0], (weighting[:,1]), axis=0, kind="linear")

    # mapping points
    numb_resp = 101
    mapping_point = np.linspace (0, 1, numb_resp)
    weight_val = smt_w (mapping_point)
    mapping_cum = np.linspace (0, 1, numb_resp)
    for ii in range (1,numb_resp) :
        mapping_cum[ii] = 0.5 * (weight_val[ii] + weight_val[ii-1])    
    mapping_cum = np.cumsum (mapping_cum)
    mapping_cum = mapping_cum / mapping_cum[-1]
    mapping = interp1d (mapping_point, mapping_cum, axis=0, kind="cubic")
    global_mapping = mapping
    
    # build new string
    alpha_eq    = np.linspace (0, 1, new_numb_node)
    alpha_new   = np.linspace (0, 1, new_numb_node)
    for ii in range(alpha_eq.shape[0]) :
        if (ii == 0) :
            alpha_new[ii] = 0
        else :
            if (ii == alpha_eq.shape[0] - 1): 
                alpha_new[ii] = 1
            else :
                global_line = alpha_eq[ii]
                alpha_new[ii] = bisect (__solve_func, 0, 1)
    return smooth_str (alpha_new)    
                     
