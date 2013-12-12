import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
import os
#pauls columns
#x rho v p ?<-dontworry about this
def L(x, a, c, N):
    con = (len(a)/len(c))
    err = 0
    dx = 1./float(N)
    for i in range(len(c)-1):
        diff = ((c[i] - a[int((i)*con)]))
        err += np.abs(diff*dx)
        
    return(err)
def makeprim(rho, mom, E):
    gam = 1.4
    eps = E - mom**2 / (2*rho)
    P = eps*(gam - 1)
    V = mom/rho
    data = [P, V]
    return(data)
Ns, L1 = [], []
filename = 'final.dat'
anal = 'solution.dat'
data = np.loadtxt(anal)
xa = data[:,0]
va = data[:,2]
rhoa = data[:,1]
Pa = data[:,3]
Nx = 5
count = 1
plt.figure()
figrho = plt.subplot(131)
plt.ylabel(r'$\rho$')
figP = plt.subplot(132)
plt.ylabel(r'Pressure')
plt.xlabel('Position')
figV = plt.subplot(133)
plt.ylabel(r'Velocity')
while(count < 8):
    Nx *= 2
    cmnd = './prob1 ' + str(Nx)
    os.system(cmnd)
    data = np.loadtxt(filename)
    lab = 'Nx= ' + str(Nx)
    x = data[:,0]
    rho = data[:,1]
    mom = data[:,2]
    E = data[:,3]
    [P, V] = makeprim(rho, mom, E)
    #print(L(x, Pa, P, Nx))
    L1.append(L(x, rhoa, rho, Nx))
    Ns.append(float(Nx))
    figrho.plot(x,rho, label = lab)
    figP.plot(x, P)
    figV.plot(x, V, label = lab)
    count +=1

figrho.plot(xa,rhoa,'--' ,label = 'exact solution')
figP.plot(xa, Pa, '--')
figV.plot(xa, va, '--')


figV.legend(loc = 'center', bbox_to_anchor = (0.5, 0.3))
'''
fit1, Ns1 = [], []
fit2, Ns2 = [], []
fit3, Ns3 = [], []
N1 = Ns[0:6]
L11 = L1[0:6]
N2 = Ns[6:11]
N3 = Ns[11:19]
L12 = L1[6:11]
L13 = L1[11:19] 


m1, b1 = np.polyfit(np.log10(N1), np.log10(L11), 1)
m2, b2 = np.polyfit(np.log10(N2), np.log10(L12), 1)
m3, b3 = np.polyfit(np.log10(N3), np.log10(L13), 1)
print(m1,b1)
print( m2, b2)
print(m3, b3)
for n in N1:
    fit1.append((10**b1 *n**m1))
    Ns1.append(n)
for n in N2:
    fit2.append((10**b2 *n**m2))
    Ns2.append(n)
for n in N3:
    fit3.append((10**b3 *n**m3))
    Ns3.append(n)

lab1 = 'slope = '+str(m1)
lab2 = 'slope = '+str(m2)
lab3 = 'slope = '+str(m3)
plt.figure()
plt.plot(Ns, L1, 'o')
plt.plot(Ns1, fit1, label = lab1)
plt.plot(Ns2, fit2, label = lab2)
plt.plot(Ns3, fit3, label = lab3)
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'L$_1$ Error')
plt.xlabel('Resolution')
plt.legend()
'''
plt.show()
