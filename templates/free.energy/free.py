#!/usr/bin/env python3

import os
import argparse
import numpy as np
import subprocess as sp
import logging
import BasisGaussian as bg
import FreeEnergyFit as fef

def func (x) :
    """ (x_1^2 -1)^2 + (x_2 + x_1^2 -1)^2 """
    tmp0 = x[0] * x[0] - 1
    tmp1 = x[1] + tmp0
    return tmp0 * tmp0 + tmp1 * tmp1        

def test_compute_force (xx):
    yy = np.zeros (np.size(xx))
    yy[0] = 8*xx[0]*(1-xx[0]*xx[0]-xx[1]*0.5)
    yy[1] = 2*(1-xx[0]*xx[0]-xx[1])
    return yy

if __name__ == "__main__" :
    bg = bg.BasisGaussian ([1.2,0], 3)
    xx = np.array ([[2., 0.], [-1.,0.]])
    hh = 1e-4
    hx = np.array ([[hh, 0], [hh,0] ])    
    # print ("value is %s" % bg.value (xx))
    # print ("value is %s" % bg.value (xx+hx))
    # print ("value is %s" % bg.value (xx-hx))
    # print ("grad  is %s" % bg.grad (xx))
    # print ("cdiff gradx is %s" % ( (bg.value (xx+hx) - bg.value(xx-hx)) / (2*hh) ))
    # print ("acc is", bg.grad (xx)[0][0] - ( (bg.value (xx+hx) - bg.value(xx-hx)) / (2*hh) )[0]) 
    # print ("acc is", bg.grad (xx)[1][0] - ( (bg.value (xx+hx) - bg.value(xx-hx)) / (2*hh) )[1]) 
    
    fe = fef.FreeEnergyFit ()
    
    nx = 30
    ny = 30
    xlo = -1.5
    xup = 1.5
    ylo = -1.5
    yup = 1.5
    hx = (xup - xlo) / nx
    hy = (yup - ylo) / ny
    xx = np.arange (xlo, xup+hx, hx)
    yy = np.arange (ylo, yup+hy, hy)    
    z_array = []
    f_array = []
    for ii in range (nx+1) :
        for jj in range (ny+1) :
            zz = np.array([xx[ii], yy[jj]])
            zf = test_compute_force (zz)
            z_array.append (zz)
            f_array.append (zf)

    fe.compute_a (1, z_array, f_array)

    f1_array = fe.value (z_array)

    fp = open ("results.out", "w")
    z0 = np.array([xx[0], yy[0]])
    f0 = func (z0)
    for ii in range (nx+1) :
        for jj in range (ny+1) :
            zz = np.array([xx[ii], yy[jj]])
            # fp.write ("%f\t%f\t  %f\t%f\n" % (xx[ii], yy[jj], func(zz) - f0, fe.value(zz) - f1))
            fp.write ("%f\t%f\t  %f\t%f\n" % (xx[ii], yy[jj], func(zz) - f0, f1_array[ii*(nx+1)+jj] - f1_array[0]))
        fp.write ("\n")

    fp.close()
