import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
import os

def plot(leg, fo, a):
    time, x, v =[], [], []
    
    for line in fo:
        t, y, w = line.split()
        time.append(t)
        x.append(y)
        v.append(w)
    plt.plot(time, x, label = leg )
    #if (a == 9):
    #   plt.plot(time, v, label = 'true')

def tits(meth):
    tit = 'error'
    if (meth == 1):
        tit = 'Forward Euler'
    elif (meth == 2):
        tit = 'Mid-Point'
    elif (meth == 3):
        tit = 'Backward Euler'
    elif (meth == 4):
        tit = 'Crank-Nichelson'
    return(tit)
        
filename = 'oscdata.dat'
n = 5.
T = 5.
meth = 2

for i in range(20):
    meth = 1
    cmnd = './osc1 ' + str(T) + ' ' + str(n) + ' ' + str(meth)
    os.system(cmnd)
    meth = 2
    cmnd = './osc2 ' + str(T) + ' ' + str(n) + ' ' + str(meth)
    os.system(cmnd)
    meth = 3
    cmnd = './osc3 ' + str(T) + ' ' + str(n) + ' ' + str(meth)
    os.system(cmnd)
    meth = 4
    cmnd = './osc4 ' + str(T) + ' ' + str(n) + ' ' + str(meth)
    os.system(cmnd)
    n*=2
    '''
    temp = "%1.1e" % (T/n)
    leg = 'h = ' + temp
    fo = open(filename)
    plot(leg, fo, i)
    n += 50
    fo.close()
meth = 1
for i in range(2):
    cmnd = './a.out ' + str(T) + ' ' + str(n) + ' ' + str(meth)
    os.system(cmnd)
    temp = "%1.1e" % (T/n)
    leg = tits(meth)
    fo = open(filename)
    #plot(leg, fo)
    meth += 1
    fo.close()

plt.ylim(-1.1, 1.1)
plt.legend(loc = 'center left', bbox_to_anchor = (1.01, 0.5))
plt.subplots_adjust(left=0.09, right=0.75)
plt.ylabel('Amplitude')
plt.xlabel('time')
plt.title(tits(meth))
plt.show()
'''

