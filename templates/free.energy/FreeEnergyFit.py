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
        t0 = time.time()
        gradMat = np.zeros ([ndata*ndata, dim])
        for ii in range (ndata) :
            for jj in range (ii+1) :
                gradMat[ii*ndata+jj] = self.bases[ii].grad (zz[jj])
        for ii in range (ndata) :
            for jj in range (ii,ndata) :
                gradMat[ii*ndata+jj] = - gradMat[jj*ndata+ii]
        t1 = time.time()
        gradMat = np.transpose (gradMat)
        t2 = time.time()
        print ("# after  gradMat  multiply: %f s trans %s s" % ((t1-t0), (t2-t1)), flush = True)        
        # assemble the matrix B
        print ("# before B", flush = True)
        BB = np.zeros ([ndata, ndata])
        CC = np.zeros (ndata)
        for kk in range (dim) :
            BB = BB + np.dot (gradMat[kk].reshape(ndata,ndata), gradMat[kk].reshape(ndata,ndata))
        # gradMat = np.zeros ([dim, ndata, ndata])
        # for ii in range (ndata) :
        #     for jj in range (ndata) :
        #         grad = self.bases[jj].grad (zz[ii])
        #         for kk in range (dim) :
        #             gradMat[kk][ii][jj] = grad[kk]
        # for kk in range (dim) :
        #     BB = BB + np.dot (gradMat[kk], gradMat[kk])
        print ("# after  B", flush = True)
        print ("# before C", flush = True)
        zft = zf.T
        for kk in range (dim) :
            CC = CC + np.dot (gradMat[kk].reshape(ndata,ndata), zft[kk])
        # for kk in range (dim) :
        #     CC = CC - np.dot (gradMat[kk], zft[kk])
        print ("# after  C", flush = True)
        # for ii in range (ndata) :
        #     for jj in range (ndata) :
        #         for kk in range (ndata) :
        #             vec0 = self.bases[kk].grad (zz[ii])
        #             vec1 = self.bases[jj].grad (zz[kk])
        #             BB[ii][jj] = BB[ii][jj] + np.dot (vec0, vec1)
        # for ii in range (ndata) :
        #     for jj in range (ndata) :
        #         CC[ii] = CC[ii] - np.dot (zf[jj], self.bases[jj].grad(zz[ii]))
        
        print ("# before solve", flush = True)
        self.a = np.linalg.solve (BB, CC)
        print ("# after  solve", flush = True)
        # print (BB)
        # print (CC)
        # print (self.bases)
        
        
