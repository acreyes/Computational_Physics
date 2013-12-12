import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
import os

def convert(t, T):
    N = len(t)
    dt = t[N-1]/N
    i = int(T/dt)
    return(i)
    

data = np.loadtxt('entropylong.dat')
t = data[:,0]
sumv = data[:,1]
sump = data[:,2]
T = 0.021
i = convert(t, 0.018)
print i
fitT = t[0:i]
fitv = sumv[0:i]
m, b = np.polyfit(fitT, fitv, 1)
loglist = []
for n in fitT:
    loglist.append(10**b *n**m)
    


plt.plot(fitT, loglist, label = 'slope = ' +str(m))
plt.plot(t, sumv, label = 'velocity diff')
plt.plot(t, sump, label = 'pressure diff')
plt.axvline(T, 0, 100, label = 't =' + str(T))
plt.legend(loc = 'lower right')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('time')
plt.ylabel('Difference')
plt.show()
