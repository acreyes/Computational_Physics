import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb

time = []
energy = []
pot = []
kin = []
filename = 'Energy.dat'
file = open(filename)
(title, temp) = filename.split('.')

for line in file:
    t, u, k, E = line.split()
    #t = float(t)
    #temp = float(E)
    

    time.append(t)
    energy.append(E)
    pot.append(u)
    kin.append(k)
   

file.close()

plt.plot(time,energy, 'b')
plt.plot(time,pot, 'r--')
plt.plot(time,kin, 'g-.')
#plt.plot(np.cos(time),np.sin(time))
#plt.plot(time, x)

plt.show()
