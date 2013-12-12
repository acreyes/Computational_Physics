import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
import os

def H(N):
    return(3/N)

def plot(leg, fo):
    time, A, X, Y =[], [], [], []
    
    for line in fo:
        t, a, x, y = line.split()
        time.append(t)
        A.append(a)
        X.append(x)
        Y.append(y)
    plt.plot(time, A,  'x',label = ('A' + leg))
    #plt.plot(time, X,  'o',label = ('X' + leg))
    #plt.plot(time, Y,  '*',label = ('Y' + leg))

    
N = 50
method = 1
filename = 'nonlindata.dat'
fo = open(filename)
for i in range(5):
    cmnd = './nonlin ' + str(method) + ' ' + str(H(N)) 
    os.system(cmnd)
    leg = 'h = '+ str(H(N))
    plot(leg, fo)
    N+=100*i
plt.legend()
plt.show()
'''
time, A, X, Y =[], [], [], []

for line in fo:
    t, a, x, y = line.split()
    time.append(t)
    A.append(a)
    X.append(x)
    Y.append(y)

plt.ylim(0, 200)
#plt.xlim(0, 0.1)
plt.plot(time, A,  'x',label = 'A')
plt.plot(time, X,  'o',label = 'X')
plt.plot(time, Y,  '*',label = 'Y')
#plt.xscale('log')
plt.legend()
fo.close()
plt.show()
'''
