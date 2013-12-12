import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

N, E, RK2, RK4 = [], [], [], []

filename = 'L2norms.dat'
file = open(filename)
for line in file:
    n, e, rk2, rk4 = line.split()
    print(n)
    N.append(n)
    E.append(e)
    RK2.append(rk2)
    RK4.append(rk4)
file.close()
plt.plot(N, E,'+', label = 'Forward Euler')
plt.plot(N, RK2,'o', label = 'Midpoint')
plt.plot(N, RK4,'x', label = 'RK4')
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Number of Steps')
plt.ylabel('L2 Norm')

plt.show()
