import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

time = []
x = []
y = []

filename = 'asign1.dat'
file = open(filename)
(title, temp) = filename.split('.')

for line in file:
    t, xtemp, ytemp = line.split()
    t = float(t)
    xtemp = float(xtemp)
    ytemp = float(ytemp)

    time.append(t)
    x.append(xtemp)
    y.append(ytemp)

file.close()

plt.plot(x,y)
plt.plot(np.cos(time),np.sin(time), 'r--')
#plt.plot(time, x)

plt.show()
