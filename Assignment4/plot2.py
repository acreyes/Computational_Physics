import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
import os
tit = (r'$\alpha$ = 1.3')
def L(x, a, c, N):
    con = (len(a)/len(c))
    err = 0
    dx = 1./float(N)
    for i in range(len(c)-1):
        diff = ((c[i] - a[int((i)*con)]))
        err += np.abs(diff*dx)
        
    return(err)
Ns, L1 = [], []
filename = 'final2.dat'
anal = 'solution.dat'
data = np.loadtxt(anal)
xa = data[:,0]
Pa = data[:,3]
Nt = 0.0
count = 1

plt.figure()
while(Nt < 0.03):
    Nt += 0.002
    cmnd = './prob2 ' + str(Nt)
    os.system(cmnd)
    data = np.loadtxt(filename)
    lab = 'T= ' + str(Nt)
    x = data[:,0]
    P = data[:,1]
    plt.plot(x,P, label = lab)
    count +=1

filename = 'initial2.dat'
data = np.loadtxt(filename)
x = data[:,0]
P = data[:,1]
plt.plot(x,P, label = 'T= 0')
plt.legend()
plt.title(tit)
plt.xlim(-.5, 0.25)
plt.xlabel('Position')
plt.ylabel(r'$\rho$')
'''
plt.figure()

filename = 'entropy.dat'
data = np.loadtxt(filename)
t = data[:,0]
vsum = data[:,1]
psum = data[:,2]
plt.plot(t, vsum, label = 'velocity difference')
plt.plot(t, psum, label = 'pressure difference')
plt.ylabel('difference')
plt.xlabel('time')
plt.yscale('log')
plt.xscale('log')
'''
plt.show()
