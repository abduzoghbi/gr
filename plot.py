#!/usr/bin/env python

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt('tmp.dat')


fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlim3d(-10,10)
ax.set_ylim3d(-10,10)
ax.set_zlim3d(-10,10)

ax.plot(d[:,0],d[:,1],d[:,2],'o')


u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

rh = 1.6
x = rh * np.outer(np.cos(u), np.sin(v))
y = rh * np.outer(np.sin(u), np.sin(v))
z = rh * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='k')


plt.show()





