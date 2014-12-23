#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
import matplotlib.cm as cm


if __name__ == '__main__':
    p           =   argparse.ArgumentParser(
        description="process the output of image_flux_time",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    p.add_argument("infile" , metavar="infile" , type=str, 
        help="input file, this is the output of gr::disk::image_flux_time")

    args = p.parse_args()
    
    fp      = open(args.infile);
    desc = fp.readline().split(); xcent = np.array(desc[2:],float)
    desc = fp.readline().split();
    taxis = desc[0] == '#'
    if taxis: tcent = np.array(desc[2:],float)


    if taxis:
        npix    =   len(xcent)
        ntime   =   len(tcent)
        images  =   np.loadtxt( args.infile ).reshape((ntime,npix,npix))
    else:
        images  =   np.loadtxt( args.infile )
    
    ind     = images==0; images[ind] += 1e10
    images  =   np.log10( images )
        
    
    ext     =   [xcent[0],xcent[-1],xcent[0],xcent[-1]]
    
    fig     =   plt.figure()
    if taxis:
        ims     =   []
        for it in range(len(tcent)):
            im =    plt.imshow( images[it] , cmap=cm.hot_r, origin='lower' , extent=ext , vmin=1,vmax=6)
            ims.append([im])
        ani = animation.ArtistAnimation(fig, ims, interval=100, repeat_delay=1000)
    else:
        plt.imshow( images , cmap=cm.hot_r, origin='lower' , extent=ext , vmin=1,vmax=6)
    
    plt.show()
    