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

np.random.seed(2017)
mean=np.zeros(6)
cov=np.eye(6)
size = 10000
epsilon = np.arange(0.2,0.725,0.025)
dim = 5
t = 1
print(size)
print(epsilon)
np.set_printoptions(precision = 6, threshold=np.inf)

def gendata(mean,cov,size):
    np.random.seed(2017)
    X = np.random.multivariate_normal(mean,cov,size)
    Xnorm = np.linalg.norm(X,axis=1) 
    Da = X/(Xnorm[:,np.newaxis])       
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
