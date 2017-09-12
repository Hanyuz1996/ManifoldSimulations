#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 22:40:12 2017

@author: mac
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 22:00:09 2017

@author: mac
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
my_path = os.getcwd()

# Simulation settings
np.random.seed(2017)
size = 20000
epsilon = np.arange(0.025,0.425,0.025)
dim = 3
t = 1
repeat_time = size
num_of_nbr = np.zeros(size)
np.set_printoptions(precision = 6, threshold=np.inf)

def gendata(size):
    u = np.asarray(np.random.uniform(0,2*np.pi,size))
    n = 0
    v = np.ndarray((size,))
    while(n<size):
        x = np.random.uniform(0,1,1)
        y = np.random.uniform(0,np.pi,1)
        fx = np.cos(y)*np.pi
        if(x<fx):
            v[n] = y
            n += 1
        else:
            continue
    Da = np.zeros((size,3))
    Da[:,0] = np.cos(u)*np.sin(v)
    Da[:,1] = np.sin(u)*np.sin(v)
    Da[:,2] = -np.cos(v)-np.log(np.tan(0.5*v))
    return(Da)


def findnbr(Da,epsilon):
    dnbr = dict()
    wnbr = []
    indnbr = []
    num_of_nbr = np.zeros(size)
    for p in range(size):
        dist = np.sqrt(np.sum((Da-Da[p])**2,axis=1))
        indp = np.where(dist < epsilon)[0]
        indp = np.setdiff1d(indp,np.array(p))
        num_of_nbr[p] = len(indp)
        for q in indp:
            dnbr[p,q] = dist[q]    
        dist = dist[indp]
        indnbr.append(indp)
        wnbr.append(np.exp(-(dist**2)/np.sqrt(epsilon)))
    return(dnbr,wnbr,indnbr,num_of_nbr)

Da = gendata(size)
meannbr = np.zeros(len(epsilon))
'''
fig = plt.figure()
ax = plt.subplot(111,projection='3d')
ax.scatter(Da[:,0],Da[:,1],Da[:,2],s=np.repeat(1,size))

plt.show()

'''
i = 0
for eps in epsilon:
    print(eps)
    temp = findnbr(Da,eps)
    meannbr[i] = np.mean(temp[3])
    i += 1
    
plt.plot(np.log10(epsilon),np.log10(meannbr))

