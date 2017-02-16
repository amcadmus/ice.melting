#!/usr/bin/env python3

import numpy as np

def string_arc_seg (string) :
    ### compute arc (alpha) segments
    numb_node = string.shape[0]
    alpha_seg = np.zeros (numb_node)
    alpha_seg[0] = 0
    for jj in range (1, numb_node):
        alpha_seg[jj] = np.linalg.norm (string[jj] - string[jj-1])
    return alpha_seg

def string_arc (string) :
    ### compute arc parameter
    alpha_seg = string_arc_seg (string)
    alpha = np.cumsum (alpha_seg)
    return alpha

def string_arc_norm (string) :
    ### compute normalized arc parameter
    alpha = string_arc (string)
    alpha = alpha / alpha[-1]
    return alpha
