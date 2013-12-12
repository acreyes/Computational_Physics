import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import pdb
import os

filename = 'relax.dat'
xx, yy, zz = np.loadtxt(filename).T
#xx = data[:,0]
#yy = data[:,1]
#zz = data[:,2]
N = np.sqrt(len(yy));
#plt.scatter(xx,yy,c=zz, s=500, marker='s')
x = np.reshape(yy, (N,N))
y = np.reshape(xx, (N,N))
z = np.reshape(zz, (N,N))
#plt.colorbar()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, z)
plt.show()
