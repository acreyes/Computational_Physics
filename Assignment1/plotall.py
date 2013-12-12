import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
filename = 'assign1_N0050_.dat\n'
time = []
X4 = []
Xe = []
X2 = []
'''
vx4 = []
vxe = []
vx2 = []
'''
time = []
Y4 = []
Ye = []
Y2 = []
'''
vy4 = []
vye = []
vy2 = []
'''
E4 = []
L4 = []

E2 = []
L2 = []

Ee = []
Le = []


trash, N, temp = filename.split('_')
#N = N.strip('N')



file = open(filename)
for line in file:
    t, xe, ye, v_xe, v_ye, x2, y2, v_x2, v_y2, x4, y4, v_x4, v_y4,Edife, Edif2, Edif4, Ldife, Ldif2, Ldif4 = line.split()

    time.append(float(t))
    Xe.append(xe)
    Ye.append(ye)
    #vxe.append(v_xe)
    #vye.append(v_ye)
    X2.append(x2)
    Y2.append(y2)
    #vx2.append(v_x2)
    #vy2.append(v_y2)
    X4.append(x4)
    Y4.append(y4)
    #vx4.append(v_x4)
    #vy4.append(v_y4)
    Ee.append(Edife)
    E2.append(Edif2)
    E4.append(Edif4)
    Le.append(Ldife)
    L2.append(Ldif2)
    L4.append(Ldif4)
file.close()

plt.figure()
plt.plot(Xe, Ye, 'r+', label = 'Forward Euler')
plt.plot(X2, Y2, 'gx', label = 'Midpoint')
plt.plot(X4, Y4, 'bo', label = 'RK4')
plt.plot(np.cos(time), np.sin(time), label = 'rtrue')
plt.xlabel('u')
plt.ylabel('w')
plt.legend()
'''
plt.figure()
plt.plot(time, Le, 'r+',label = 'Forward Euler')
plt.plot(time, L2, 'gx', label = 'Midpoint')
plt.plot(time, L4, 'bo', label = 'RK4')
plt.xlabel('time')
plt.ylabel('Angular Momentum')
plt.legend()

plt.figure()
plt.plot(time, Ee,label = 'Forward Euler')
plt.plot(time, E2, label = 'Midpoint')
plt.plot(time, E4, label = 'RK4')
plt.yscale('log')
plt.xlabel('time')
plt.ylabel('Total Energy')
plt.legend()
'''
plt.show()





