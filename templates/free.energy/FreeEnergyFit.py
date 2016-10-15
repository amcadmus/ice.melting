#!/usr/bin/env python3

import numpy as np
import BasisGaussian as bg
import time

class FreeEnergyFit (object) :
    """
    Fit the free energy profile based on forces
    """
    def __init__ (self) :
        self.bases = []
        self.a = []
        return
    
    def value (self,
               xx ):
        value = 0
        for ii in range(len(self.bases)) :
            value = value + self.a[ii] * self.bases[ii].value(xx)
        return value

    def compute_a (self,
                   sigma,
                   zz_,
                   zf_ ) :
        """
        compute the prefactors a given a sigma
        zz: the position
        zf: the force measured at zz
        """
        zz = np.array(zz_)
        zf = np.array(zf_)
        sigma = zz[1][1] - zz[0][1] 
        ndata = zz.shape[0]
        dim = zz.shape[1]
        self.bases = []
        for ii in range (ndata) :
            self.bases.append (bg.BasisGaussian(zz[ii], sigma))

        # generate the matrix gradMat
        print ("# before gradMat", flush = True)
        tm1 = time.time()
        gradMat = np.zeros ([ndata*ndata, dim])
        for ii in range (ndata) :
            gradMat[ii*ndata:ii*ndata+ii+1] = self.bases[ii].grad (zz[0:ii+1]) 
            # for jj in range (ii+1) :
            #     gradMat[ii*ndata+jj] = self.bases[ii].grad (zz[jj]) 
            # tmp = self.bases[ii].grad (zz[0:ii+1]) 
            # if (np.linalg.norm (tmp - gradMat[ii*ndata:ii*ndata+ii+1])) > 1e-12 : 
            #     print (tmp)
            # print (tmp)
            # print (gradMat[ii*ndata:ii*ndata+ii+1])
        t0 = time.time()
        for ii in range (ndata) :
            for jj in range (ii,ndata) :
                gradMat[ii*ndata+jj] = - gradMat[jj*ndata+ii]
        t1 = time.time()
        gradMat = np.transpose (gradMat)
        gradMat = gradMat.reshape (dim, ndata, ndata)
        t2 = time.time()
        print ("# after  gradMat  total %f s, multiply0: %f s, multiply1: %f s, trans %s s" % ( (t2 - tm1), (t0 - tm1), (t1-t0), (t2-t1)), flush = True)        

        # assemble the matrix B
        BB = np.zeros ([ndata, ndata])
        CC = np.zeros (ndata)
        print ("# before B", flush = True)
        t0 = time.time()
        for kk in range (dim) :
            BB = BB + np.matmul (gradMat[kk], gradMat[kk])
        t1 = time.time()
        print ("# after  B   total %f s" % (t1 - t0), flush = True)

        # assemble the matrix C
        print ("# before C", flush = True)
        zft = zf.T
        for kk in range (dim) :
            CC = CC + np.dot (gradMat[kk], zft[kk])
        print ("# after  C", flush = True)
        
        # assemble the matrix C
        print ("# before solve", flush = True)
        t0 = time.time()
        self.a = np.linalg.solve (BB, CC)
        t1 = time.time()
        print ("# after  solve   total %f s" % (t1 - t0), flush = True)
        # print (BB)
        # print (CC)
        # print (self.bases)
        
        
