import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

filename = 'assign2.dat'
fo = open(filename)
time, x4, y4, xv, yv, E4, Ev = [], [], [], [], [], [], []

for line in fo:
    t, x1, y1, x2, y2, e4, ev = line.split()
    time.append(t)
    x4.append(x1)
    y4.append(y1)
    xv.append(x2)
    yv.append(y2)
    E4.append(e4)
    Ev.append(ev)

#plt.subplot(211)
#plt.plot(x4, y4, 'ro', label = 'RK4')
plt.plot(xv, yv, 'b-', label = 'Verlet')
plt.legend()
#plt.text(-3.5,1, 'N = 20,000\nT=6000')
#plt.xlim(-4,2)
#plt.ylim(-4,2)
plt.xlabel('x')
plt.ylabel('y')

'''
plt.subplot(212)
plt.plot(time, E4, 'ro')
plt.plot(time, Ev )
plt.xlabel('time')
plt.ylabel(r'$|E - E_0|/E_0$')
plt.yscale('log')
'''
plt.show()
