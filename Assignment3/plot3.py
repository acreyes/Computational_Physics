import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

def f(t):
    return(1-10*t**2)

filename = 'oscdata.dat'
fo = open(filename)
time, x, v =[], [], []

for line in fo:
    t, y, w = line.split()
    time.append(t)
    x.append(y)
    v.append(w)

plt.plot(time, x, 'r', label = 'high res Mid-Point')
plt.plot(time, v, label = 'Analytic Solution')
plt.xlabel('time')
plt.ylabel('Amplitude')
fo.close()
plt.legend()
plt.show()
