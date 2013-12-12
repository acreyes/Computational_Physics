import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

filename = 'l2.dat'
fid = open(filename)
N, L, fit = [], [], []
def f(x, m, b):
    return(m*x + b)
for line in fid:
    n, l = line.split()
    N.append(np.log10(float(n)))
    L.append(np.log10(float(l)))

m, b = np.polyfit(N, L, 1)
print(m,b)
for n in N:
    fit.append(m*n + b)
plt.plot(N, L, 'ro')
plt.plot(N, fit)
plt.text(1.5, -2, 'slope = -0.53')
plt.ylabel(r'$log_{10} L2$')
plt.xlabel(r'$log_{10} N$')
#plt.yscale('log')
#plt.xscale('log')
plt.show()
