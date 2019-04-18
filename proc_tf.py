#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import filters


if __name__ == '__main__':
    p           =   argparse.ArgumentParser(
        description="process the output of gr::disk::tf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    p.add_argument("infile" , metavar="infile" , type=str, 
        help="input file, this is the output of gr::disk::tf")

    args = p.parse_args()
    
    fp      = open(args.infile);
    desc    = fp.readline().split(); tcent = np.array(desc[2:],float)
    desc    = fp.readline().split(); ecent = np.array(desc[2:],float)
    
    
    tf      = np.loadtxt( args.infile )
    idx = tcent<80
    tf = tf[:, idx]
    tcent = tcent[idx]
    #tf      = filters.gaussian_filter(tf, 0.5)
    
    ind     = tf==0
    #tf[ind] += 1e10
    #tf      = np.log10(tf)

        
    
    ext     =   [tcent[0],tcent[-1],ecent[0],ecent[-1]]
    
    
    plt.imshow(np.log(tf) , cmap=cm.hot_r, origin='lower' , extent=ext  , aspect='auto')
    #plt.show();exit(0)
    from IPython import embed;embed();exit(0)
    plt.subplot(121); plt.plot(ecent, tf.sum(1))
    plt.subplot(122); plt.plot(tcent, tf.sum(0))
    plt.show()

    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')

    
