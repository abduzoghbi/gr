#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import matplotlib.cm as cm


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
    
    ind     = tf==0; tf[ind] += 1e10
    tf      = np.log10(tf)

        
    
    ext     =   [tcent[0],tcent[-1],ecent[0],ecent[-1]]
    
    fig     =   plt.figure()
    plt.imshow( tf , cmap=cm.hot_r, origin='lower' , extent=ext  , aspect='auto')
    
    plt.show()
    