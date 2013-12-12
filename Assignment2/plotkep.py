import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

filename = 'kepler.dat'
name = 'assign2.dat'
fid = open(name)
fo = open(filename)

time, x, y, xv, yv = [], [], [], [], []

for line in fid:
    t, x1, y1, x2, y2, e4, ev = line.split()
    xv.append(x2)
    yv.append(y2)

for line in fo:
    t, a, b = line.split()
    time.append(t)
    x.append(a)
    y.append(b)


plt.figure()
plt.plot(x, y, label = 'Kepler')
plt.plot(xv, yv, label = 'Verlet')
plt.xlabel('x')
plt.ylabel('y')
plt.text(-2, 1., 'T = 200\nN = 5000')
plt.legend()
plt.show()
