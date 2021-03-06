# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
O"""
import time
import numpy as np
import numpy.linalg as npl
import scipy as sp
import scipy.stats as st
from scipy.spatial.distance import pdist,squareform
#import matplotlib.pyplot as plt
#import pandas as pd
import math
import os
my_path = os.getcwd()

# Simulation settings
np.random.seed(2017)
mean=[0,0,0]
cov=[[1,0,0],[0,1,0],[0,0,1]]
size = 2000
epsilon = np.array([0.3])
d = 2
t = 1
sample_time = size/10
repeat_time = size
num_of_nbr = np.zeros(size)
print(size)
print(epsilon)

# Generating uniformly distributed data on a 2D sphere
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
    for p in range(size):
        dist = np.sqrt(np.sum((Da-Da[p])**2,axis=1))
        indp = np.where(dist < epsilon)[0]
        indp = np.setdiff1d(indp,np.array(p))
        for q in indp:
            dnbr[p,q] = dist[q]
        dist = dist[indp]
        indnbr.append(indp)
        wnbr.append(np.exp(-(dist**2)/epsilon))
    return(dnbr,wnbr,indnbr)

def findbase(Da,dnbr,wnbr,indnbr):
    Oplist = []
    for p in range(size):
        if len(indnbr[p]) < 2:
            Oplist.append(np.array([-1]))
            continue
        Xp = (Da[indnbr[p]] - Da[p]).T # Only near neighbor is used 
        Dp = np.diag(np.sqrt(wnbr[p]))    
        Bp = np.dot(Xp, Dp)
        u, sigma, v = npl.svd(Bp)
        Op = u[:,:2]
        Oplist.append(Op)
    '''
    Q = dict()
    for i in range(size):
        for j in range(size):
            Q[i,j]=np.dot(Oplist[i].T, Oplist[j])
    B = np.zeros((size*d,size*d))
    for i in range(size):
        for j in range(size):
            temp = np.zeros((2,2))
            if j==i:
                for k in indnbr[p]:
                    temp = (np.eye(d) + np.dot(Q[j,k],Q[j,k].T))*wmat[j,k]
                    B[i*d:(i+1)*d,:][:,j*d:(j+1)*d] = B[i*d:(i+1)*d,:][:,j*d:(j+1)*d] + temp
            else:
                B[i*d:(i+1)*d,:][:,j*d:(j+1)*d] = -2*wmat[i,j]*Q[i,j]
    ''' 
    return(Oplist)
  
def geodir(p,Da,Oplist,indnbr,d):
    Base = Oplist[p].T
    Xp = (Da[indnbr[p]] - Da[p]).T
    v0 = np.dot(Base, Xp)
    vnorm = npl.norm(v0,axis=0)
    v0 = v0/vnorm.T
    Vinit = np.zeros((d,size))
    Vinit[:,indnbr[p]]=v0
#    Vinit = np.reshape(Vinit,(size*d,1),order = 'F')
#    V = npl.solve((np.eye(size*d)+t*B),Vinit)
#    V = np.reshape(V,(d,size),order = 'F')
#    vnorm = npl.norm(V,axis=0)
#    V = V/vnorm.T
    return(Vinit)

# Function findopq aims to find the parallel transport operator Opq for each 
# pair of points p,q in the data point. The input is the base at each point 
# Oplist and epsilon. To accelarate the computation, we only compute this 
# parallel transport operator for near neighbors in the epsilon ball
'''
def findopq(Oplist, epsilon):
    Opq = dict()
    for p in range(size):
        for q in range(size):
            if dmat[p,q]>=epsilon:
                continue
            else:
                Op = Oplist[p]
                Oq = Oplist[q]
                u,_,v = npl.svd(np.dot(Op.T, Oq))
                Opq[p, q] = np.dot(u, v)
    return(Opq)
'''

def findopq(Oplist,epsilon,dnbr):
    Opq = dict()
    for p in range(size):
        if (Oplist[p] == -1).all():
            continue
        for q in range(size):
            if (Oplist[q] == -1).all():
                continue
            elif (p,q) not in dnbr:
                continue
            else:
                Op = Oplist[p]
                Oq = Oplist[q]
                u,_,v = npl.svd(np.dot(Op.T, Oq))
                Opq[p, q] = np.dot(u, v)
    return(Opq)

# Function curvatp aims to compute the Gaussian curvature in this sphere case.
# We already obtain the orthonormal basis for tangent plane at every point.
# Therefore in this local base the coordinate of the tangent vector is [1,0] and
# [0,1]. We compute R(X,Y,Y,X) with X=[1,0],Y=[0,1] by first computing R(X,Y)Y
# in our approach of parallel transport and triangle approximation. Then we obtain
# the gaussian curvature by compute the inner product of R(X,Y)Y and X, which means
# we only need to select the first item in the result of R(X,Y)Y. This function
# only gives the result of R(X,Y)Y.

# A significant problem I didn't notice before is the signal of the curvature
# Here to test the effect of the estiamtion, I naively save the absolute value
# of the estimation result since we know that the Gaussian curvature for a sphere
# should be positive.
def curvatp_paratrans_simp(p,Da,Oplist,Opq,minind,dnbr):
    V = geodir(p,Da,Oplist,indnbr,d)
    if((minind[1],minind[2]) not in dnbr):
        return(-1)
    else:
        u = np.random.uniform(low=0,high=2*np.pi,size=1)
        X=-V[:,minind[1]]
        Y=-V[:,minind[2]]
        pz = np.dot(Opq[minind[0],minind[1]],np.dot(Opq[minind[1],minind[2]],
                    np.dot(Opq[minind[2],minind[0]],Y)))
        l1 = dnbr[minind[0],minind[1]]
        l2 = dnbr[minind[1],minind[2]]
        l3 = dnbr[minind[0],minind[2]]
        s = 0.5 * (l1 + l2 + l3)
        Area = math.sqrt(s*(s-l1)*(s-l2)*(s-l3))
        diff = ((pz - Y)/Area)
        res = np.inner(diff,X)/(1-np.inner(X,Y)*np.inner(X,Y))
#        print("R tensor: ", np.inner(diff,X), "Area: ", Area)
    return(abs(res)) 
    
def curvatp_paratrans_ensemble(p,Da,Oplist,Opq,epsilon,dnbr,indnbr):
    nbr = indnbr[p]
    if(len(nbr)<2):
        res = -1
        return(res)
    K = np.min([sample_time,len(nbr)*(len(nbr)-1)/2])
    sample_res = np.zeros(K)
    for time in range(K):
        ind = np.random.choice(nbr,2,replace = False)
        ind = np.append(p,ind)
        sample_res[time] = curvatp_paratrans_simp(p,Da,Oplist,Opq,ind,dnbr)
    if(len(sample_res[sample_res>0])==0):
        return(-1)
    else:
        res = np.median(sample_res[sample_res>0])
    return(res)
    
def angSum(i,j,k,Oplist,Da):
    basei = Oplist[i].T
    basej = Oplist[j].T
    basek = Oplist[k].T
    iij = np.dot(basei,(Da[j]-Da[i]).T)
    iij = iij/npl.norm(iij)
    iik = np.dot(basei,(Da[k]-Da[i]).T)
    iik = iik/npl.norm(iik)
    jji = np.dot(basej,(Da[i]-Da[j]).T)
    jji = jji/npl.norm(jji)
    jjk = np.dot(basej,(Da[k]-Da[j]).T)
    jjk = jjk/npl.norm(jjk)
    kki = np.dot(basek,(Da[i]-Da[k]).T)
    kki = kki/npl.norm(kki)
    kkj = np.dot(basek,(Da[j]-Da[k]).T)
    kkj = kkj/npl.norm(kkj)
    angi = np.arccos(np.inner(iij,iik))
    angj = np.arccos(np.inner(jji,jjk))
    angk = np.arccos(np.inner(kki,kkj))
    return(angi+angj+angk)
    
def findArea(i,j,k,dnbr):
    if(((i,j) not in dnbr) or ((i,j) not in dnbr) or ((k,j) not in dnbr)):
        return(-1)
    else:
        l1 = dnbr[i,j]
        l2 = dnbr[j,k]
        l3 = dnbr[k,i]
        s = 0.5*(l1+l2+l3)
        area = np.sqrt(s*(s-l1)*(s-l2)*(s-l3))
    return(area)
 
def curvatp_gb_simp(p,ind,Oplist):
    if(p in ind):
        ind = np.setdiff1d(ind,np.array(p))
        ang = angSum(p,ind[0],ind[1],Oplist,Da)
        Area = findArea(p,ind[0],ind[1],dnbr)
        if(Area < 0):
            print(p, "not in neighborhood")
            return(-3)
        else: 
            res = (ang-np.pi)/Area
    else:
        ang = angSum(ind[0],ind[1],ind[2],Oplist,Da)
        Area = findArea(p,ind[0],ind[1],dnbr)+findArea(p,ind[1],ind[2],dnbr)+findArea(p,ind[0],ind[2],dnbr)
        if(Area < 0):
            return(-3)
        else: 
            res = (ang-np.pi)/Area
    return(res) 
'''        
def curvatp_quadfit(p,indnbr,Da,dmat,wmat,d, base):
    nbr = indnbr[p]
    X = (Da[nbr] - Da[p]).T
    norm = np.cross(base[:,0],base[:,1])
    base = np.vstack((base[:,0],base[:,1],norm))
    X = np.dot(base.T, X).T
    print("X: ",X)
    z= X[:,2]
    xsquare = X[:,0] * X[:,0]
    xy = X[:,0] * X[:,1]
    ysquare = X[:,1] * X[:,1]
    X = np.vstack((xsquare,xy,ysquare)).T
    print("X: ", X)
    beta = npl.solve(np.dot(X.T,X),np.dot(X.T,z))
    print("beta: ", beta)
    res = 4*beta[0]*beta[2]-beta[1]*beta[1]
    return(res)    
'''

def curvatp_gb_enclose(p,indnbr,Da,Oplist,dnbr):
    nbr = indnbr[p]
    if(len(nbr) < 3):
        return(-1)
    used = np.zeros(len(nbr))
    base = Oplist[p].T
    X = np.dot(base,(Da[nbr]-Da[p]).T)
    res_all = []
    ind = 0
    for i in np.arange(ind,len(nbr)-2):
        if used[i] == 1:
            continue
        for j in np.arange(ind+1,len(nbr)-1):
            if used[j] == 1:
                continue
            for k in np.arange(j + 1,len(nbr)):
                if used[k] == 1:
                    continue
                else:
                    norm0 = np.cross(X[:,i],X[:,j])
                    norm1 = np.cross(X[:,j],X[:,k])
                    norm2 = np.cross(X[:,k],X[:,i])
                    if(np.inner(norm0,norm1) < 0 or np.inner(norm1,norm2) < 0):
                        continue
                    else:
                        index = np.array([nbr[i],nbr[j],nbr[k]])
                        result = curvatp_gb_simp(p,index,Oplist)
                        res_all.append(result)
                        used[i] = 1
                        used[j] = 1
                        used[k] = 1
                        ind += 1
        continue
    res_all = np.array(res_all)
    res_all = res_all[res_all>0]
    if(len(res_all) == 0):
        return(-2)
    else:
        res = np.median(res_all)
    return(res)
    
#reseps_verysmall = np.zeros(50)
#reseps_25 = np.zeros(50)
#reseps_75 = np.zeros(75)
#reseps_verylarge = np.zeros(50)  
reseps_med = np.zeros(50)
reseps_mean = np.zeros(50)
reseps_var = np.zeros(50)
res_gb = np.zeros((len(epsilon),size))
Da = gendata(mean,cov,size)
iter = 0

for eps in epsilon:
    dnbr,wnbr,indnbr = findnbr(Da,eps)
    Oplist = findbase(Da,dnbr,wnbr,indnbr)
    Opq = findopq(Oplist,eps,dnbr)
    for p in range(repeat_time):
        res_gb[iter,p] = curvatp_gb_enclose(p,indnbr,Da,Oplist,dnbr)
#       res_paras_ensemble[iter,p] = curvatp_paratrans_ensemble(p,Da,Oplist,Opq,eps,dnbr,indnbr)
#       res_quadfit[p] = curvatp_quadfit(p,indnbr,Da,dmat,wmat,d, Oplist[p])
    iter += 1
print(res_gb)
   
 



    
