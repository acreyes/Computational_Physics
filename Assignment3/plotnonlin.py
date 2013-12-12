import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
import os
n = 5
meth = 1
if(meth == 1):
    title = 'Forward Euler'
else:
    title = 'Backward Euler'

filename = 'beun.dat'
file2 = 'best.dat'
fo = open(filename)
data = np.loadtxt(file2)
time2 = data[:,0]
A2 = data[:,1]
X2 = data[:,2]
Y2 = data[:,3]
time, A, X, Y =[], [], [], []
h1 = '5.e-4'
h2 = '5.e-5'
for line in fo:
    t, a, x, y = line.split()
    time.append(t)
    A.append(a)
    X.append(x)
    Y.append(y)
        
plt.subplot(122)
plt.plot(time, A,  'x',label = ('A'),ms = 2)
plt.plot(time, X,  'go',label = ('X'),ms = 2)
plt.plot(time, Y,  'r+',label = ('Y'),ms = 4)
plt.ylim(0, 155)
plt.title('Backward Euler')
plt.xscale('log')

plt.ylabel('Amount of Species')
plt.xlabel('time')
plt.legend()
plt.subplot(121)
plt.plot(time2, A2,  'x',label = ('A'),ms = 2)
plt.plot(time2, X2,  'go',label = ('X'),ms = 2)
plt.plot(time2, Y2,  'r+',label = ('Y'),ms = 2)
plt.title('Forward Euler')



plt.ylim(0, 155)
#plt.xlim(0, 0.1)
plt.xscale('log')
plt.xlabel('time')
#plt.title(r'Forward Euler with step size $10^{-5}$')
plt.legend()
fo.close()
plt.show()

