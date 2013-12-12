import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

file1 = 'L21.dat'
file2 = 'L22.dat'
file3 = 'L23.dat'
file4 = 'L24.dat'
#data = np.loadtxt('oscdata.dat')
#t = data[:,0]
#f = data[:,2]
data = np.loadtxt(file1)
h1 = data[:,0]
L1 = data[:,1]
data = np.loadtxt(file2)
h2 = data[:,0]
L2 = data[:,1]
data =  np.loadtxt(file3)
h3 = data[:,0]
L3 = data[:,1]
data = np.loadtxt(file4)
h4 = data[:,0]
L4 = data[:,1]


plt.figure()
plt.plot(h1, L1,'--', label = 'Forward Euler')
plt.plot(h2, L2,'*', label = 'MidPoint')
plt.plot(h3, L3,'o', label = 'Backward Euler')
plt.plot(h4, L4,'x', label = 'Crank Nichelson')
plt.xscale('log')
plt.yscale('log')
plt.title(r'$L_2$ Error')
plt.ylabel(r'$L_2$ Error')
plt.xlabel('Number of Steps')
plt.legend()
#plt.plot(t, f)
plt.show()
