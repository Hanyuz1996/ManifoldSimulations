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
import matplotlib.pyplot as mlp
import os
my_path = os.getcwd()

# Simulation settings
np.random.seed(2017)
mean=[0,0]
cov=[[1,0],[0,1]]
size = 20000
epsilon = np.arange(0.025,0.725,0.025)
dim = 3
t = 1
repeat_time = size
num_of_nbr = np.zeros(size)
print(size)
print(epsilon)
np.set_printoptions(precision = 6, threshold=np.inf)

def gendata(mean,cov,size):
    np.random.seed(2017)
    X1 = np.random.multivariate_normal(mean,cov,size)
    Xnorm = np.linalg.norm(X1,axis=1) 
    X1 = X1/(Xnorm[:,np.newaxis])
    X2 = np.random.multivariate_normal(mean,cov,size)
    Xnorm = np.linalg.norm(X2,axis=1) 
    X2 = X2/(Xnorm[:,np.newaxis])
    X3 = np.random.multivariate_normal(mean,cov,size)
    Xnorm = np.linalg.norm(X3,axis=1) 
    X3 = X3/(Xnorm[:,np.newaxis])
    Da = (np.vstack([X1[:,0],X1[:,1],X2[:,0],X2[:,1],X3[:,0],X3[:,1]])).T         
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

Da = gendata(mean,cov,size)
meannbr = np.zeros(len(epsilon))
i = 0
for eps in epsilon:
    print(eps)
    temp = findnbr(Da,eps)
    meannbr[i] = np.mean(temp[3])
    i += 1
