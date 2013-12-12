import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

def h(t):
    return(t**3 - np.cos(t))
def dfdx(t):
    return(3*t**2 + np.sin(t))

filename = 'newraph.dat'
fo = open(filename)
t = np.linspace(-6, 6, 300)
N, x, f = [], [], []

for line in fo:
    n, a, g = line.split()
    N.append(n)
    x.append(float(a))
    f.append(float(g))

plt.figure()
plt.plot(N, f)
plt.title('Function Evaluated at Computed Root vs # of Iterations')
plt.xlabel('# of Iterations')
plt.ylabel(r'$x^3 - cos(x)$')
#plt.yscale('log')
#plt.xscale('log')

#plt.plot(t, t**3 - np.cos(t))
plt.show()
