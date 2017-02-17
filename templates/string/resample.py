#!/usr/bin/env python3

import numpy as np
import StringUtils

if __name__ == "__main__":
    string = np.loadtxt ("string.out")
    weight = np.loadtxt ("weight.out")
    print (string.shape[0] + 10)
    new_string = StringUtils.resample_string (string, string.shape[0] + 16, weight)
    new1_string = StringUtils.resample_string (new_string, new_string.shape[0], weight)
#    new1_string = resample_string (new_string, weight, 1 + 2 * (new_string.shape[0] - 1))
    np.savetxt ("new.out", new_string)
    np.savetxt ("new1.out", new1_string)
   
