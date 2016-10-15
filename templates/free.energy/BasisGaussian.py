#!/usr/bin/env python3

import math
import numpy as np

class BasisGaussian (object) :
    """
    The Gaussian basis e^{-|x - x_c|^2 / 2}
    """
    def __init__ (self,
                  center, 
                  sigma) :
        self.center = np.array(center)
        self.dim = np.size(center)
        assert (np.size(sigma) == 1)
        self.sigma2 = sigma * sigma

    def value (self,
               xx_) :
        xx = np.array(xx_)
        uu = (np.linalg.norm (xx - self.center))
        return math.exp (- 0.5 * uu * uu / self.sigma2)

    def grad (self,
              xx_) :
        xx = np.array(xx_)
        uu = (np.linalg.norm (xx - self.center)) / self.sigma2
        return self.value (xx) * (self.center - xx) / self.sigma2
        
    
