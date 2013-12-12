import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import pdb
import os

'''
filename = 'initial.dat'
data = np.loadtxt(filename)
x = data[:,0]
y = data[:,1]
T = data[:,2]
xx, yy = np.meshgrid(x, y)
plt.figure()
plt.scatter(x, y, c=T, marker='s', s=100)
axes().set_aspect('equal', 'datalim')
plt.colorbar()
plt.title('T = 0')
'''
fig = plt.figure()
xx, yy, TT = np.loadtxt('initial.dat').T
N = np.sqrt(len(yy))
x = np.reshape(xx, (N,N))
y = np.reshape(yy, (N,N))
T = np.reshape(TT, (N,N))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, T)
#ax.set_zlim([0,12])
#axes().set_aspect('equal', 'datalim')
plt.title('T = 0')
plt.savefig(str(0))
#plt.show()

t=0.0001
dt = 0.0002
count = 1
while(count<10):
    t+= dt
    cmnd = './heateq ' + str(t) + ' 0.0001'
    os.system(cmnd)
    filename = 'heat.dat'
    xx,yy,TT = np.loadtxt(filename).T
    N = np.sqrt(len(xx))
    x = np.reshape(xx, (N,N))
    y = np.reshape(yy, (N,N))
    T = np.reshape(TT, (N,N))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x, y, T)
    #ax.set_zlim([0,12]) 
    #plt.scatter(x,y,c=T, s=500, marker='s', vmin = 0, vmax = 10)
    #axes().set_aspect('equal', 'datalim')
    #plt.colorbar()
    plt.title('T = '+str(t))
    plt.savefig(str(count))
    count+=1
    #plt.show()



'''
20 - s=500 is good

'''
