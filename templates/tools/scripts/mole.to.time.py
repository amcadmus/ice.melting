#!/usr/bin/env python3

import os
import argparse
import numpy as np
import subprocess as sp
import glob

def main () :
    parser = argparse.ArgumentParser(
        description="*** convert mole traj files matches dir/mol_* to step info in the same directory. ***")
    parser.add_argument("INPUT",
                        help='The folder of input mol trajs.')
    parser.add_argument('-c', '--column', default = 3, type=int,
                        help='The column used in the mol trajs.')

    args = parser.parse_args()
    myin = args.INPUT
    col = args.column

    if not os.path.exists (myin) :
        raise RuntimeError ("no dir " + myin)

    files = glob.glob (myin + "/mol_*")
    files.sort()
    nmols = len(files)
    print ("# numb mole %d" % nmols)

    # get steps
    data = np.loadtxt (files[0])
    steps = data[:,0].astype(int)
    nsteps = len(steps)
    print ("# numb time step %d" % nsteps)
    
    # load data in to matrix
    matrix=np.zeros ((nmols, nsteps))
    count = 0
    for ii in files :
        print ("# load from file %s" % ii, end='\r')
        data = np.loadtxt(ii)[:,col-1]
        assert (len(data) == nsteps)
        matrix[count] = data
        count = count+1
    print ("")

    # print
    assert (nsteps == matrix.shape[1])
    mat_print = np.zeros ((nmols,2))
    mat_print[:,0] = range (0, nmols)
    for ii in range (0,nsteps) :
        filename = myin + "/step_%09d" % steps[ii]
        print ("# print to  file %s" % filename, end='\r')
        mat_print[:,1] = matrix[:,ii]
        np.savetxt (filename, mat_print)
    print ("")
    

if __name__ == '__main__' : 
    main ()
           
