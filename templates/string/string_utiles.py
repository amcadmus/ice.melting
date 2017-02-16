#!/usr/bin/env python3

import numpy as np

def string_param (string)
### compute arc (alpha)
    numb_node = string.shape[0]
    alpha       = np.zeros (numb_node)
    alpha_seg   = np.zeros (numb_node)
    alpha_seg[0] = 0
    for jj in range (1, numb_node):
        alpha_seg[jj] = np.linalg.norm (string[jj] - string[jj-1])
    alpha = np.cumsum (alpha_seg)
    alpha = alpha / alpha[-1]
    return alpha
