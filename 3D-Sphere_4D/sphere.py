import time
import numpy as np
import numpy.linalg as npl
import math
import os
import sys
my_path = os.getcwd()

# Simulation settings
np.random.seed(2017)
mean=np.zeros(5)
cov=np.eye(5)
size = 20000
epsilon = 0.4
dim = 4
t = 1
repeat_time = 1000
np.set_printoptions(precision = 6, threshold=np.inf)


# Generating the stretching matrix S
def findind(typeind,i,j=1,k=1,l=1,d=dim):
    if(typeind == 1):
        res = (i-1)*(d**3)+(j-1)*(d**2)+(k-1)*d+l
    if(typeind == 2):
        res = 0.5*(2*d-i)*(i-1)+(j-i)
    if(typeind == 3):
        res = (3*d*d*i-3*d*d-3*d*i*i+3*d+i**3-i)/6-0.5*(i-j+1)*(2*d-i-j)-j+k
    if(typeind == 4):
        res = -(i-j+1)*(3*d**2-3*d*i-3*d*j-3*d+i**2+i*j+2*i+j**2+j)/6-0.5*(j-k+1)*(2*d-j-k)-k+l
        +(4*d**3*i-4*d**3-6*d**2*i**2-6*d**2*i+12*d**24*d*i**3+6*d*i**2-2*d*i-8*d-i**4-2*i**3+i**2+2*i)
    return(int(res))
                   

def stretch(d):
    S = np.zeros((d**4,d*d*(d**2-1)/12))
    for i in range(1,d+1):
        for j in range(1,d+1):
            if(i==j):
                continue
            for k in range(1,d+1):
                for l in range(1,d+1):
                    if(k==l):
                        continue
                    if(i==k):
                        if(j==l):
                            S[findind(1,i,j,k,l)-1,findind(2,np.min([i,j]),np.max([i,j]))-1]=1
                            continue
                        else:
                            sort = np.sort([i,j,l])
                            if(i == sort[0]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)]=-1
                                continue
                            elif(i == sort[1]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+1]=-1
                                continue
                            else:
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+2]=-1
                                continue
                    if(i==l):
                        if(j==k):
                            S[findind(1,i,j,k,l)-1,findind(2,np.min([i,j]),np.max([i,j]))-1]=-1
                            continue
                        else:
                            sort = np.sort([i,j,k])
                            if(i == sort[0]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)]=1
                                continue
                            elif(i == sort[1]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+1]=1
                                continue
                            else:
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+2]=1
                                continue
                    if(j==k):
                        sort = np.sort([i,j,l])
                        if(j == sort[0]):
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)]=1
                            continue
                        elif(j == sort[1]):
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+1]=1
                            continue
                        else:
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+2]=1
                            continue
                    if(j==l):
                        sort = np.sort([i,j,k])
                        if(j == sort[0]):
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)]=-1
                            continue
                        elif(j == sort[1]):
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+1]=-1
                            continue
                        else:
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+3*(findind(3,sort[0],sort[1],sort[2])-1)+2]=-1
                            continue
                    arg = np.argsort(np.argsort([i,j,k,l]))
                    sort = np.sort([i,j,k,l])
                    if(abs(arg[0]-arg[1])==abs(arg[2]-arg[3])):
                        if(abs(arg[0]-arg[1])==1):
                            if(arg[0]<arg[1] and arg[2]<arg[3]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=1
                                continue
                            elif(arg[0]<arg[1] and arg[2]>arg[3]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=-1
                                continue
                            elif(arg[0]>arg[1] and arg[2]<arg[3]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=-1
                                continue
                            else:
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=1
                                continue
                        else:
                            if(arg[0]<arg[1] and arg[2]<arg[3]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)+1]=1
                                continue
                            elif(arg[0]<arg[1] and arg[2]>arg[3]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)+1]=-1
                                continue
                            elif(arg[0]>arg[1] and arg[2]<arg[3]):
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)+1]=-1
                                continue
                            else:
                                S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)+1]=1
                                continue
                    else:
                        if(arg[0]<arg[1] and arg[2]<arg[3]):
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=-1
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)+1]=1
                            continue
                        elif(arg[0]<arg[1] and arg[2]>arg[3]):
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=1
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)+1]=-1
                            continue
                        elif(arg[0]>arg[1] and arg[2]<arg[3]):
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=1
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)+1]=-1
                            continue
                        else:
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=-1
                            S[findind(1,i,j,k,l)-1,d*(d-1)/2+d*(d-1)*(d-2)/2+2*(findind(4,sort[0],sort[1],sort[2],sort[3])-1)]=1
                            continue
    return(S)    
                    
# Generating data points
def gendata(mean,cov,size):
    X = np.random.multivariate_normal(mean,cov,size)
    Xnorm = np.linalg.norm(X,axis=1) 
    Da = X/(Xnorm[:,np.newaxis])
    Da[0] = np.array([1,0,0,0,0])
    dist_temp = np.sqrt(np.sum((Da-Da[0])**2,axis=1))
    indp_temp = np.where(dist_temp < 2*epsilon)[0]
    Da = Da[indp_temp]       
    return(Da)

S = stretch(dim)

def findnbr(Da,epsilon):
    dnbr = dict()
    wnbr = []
    indnbr = []
    indtri = []
    for p in range(len(Da)):
        dist = np.sqrt(np.sum((Da-Da[p])**2,axis=1))
        indp = np.where(dist < epsilon)[0]
        indp = np.setdiff1d(indp,np.array(p))
        num_of_nbr[p] = len(indp)
        for q in indp:
            dnbr[p,q] = dist[q]    
        indnbr.append(np.random.choice(indp,len(indp)/2,replace=False))
        indtri.append(np.setdiff1d(indp,indnbr[p]))
        wnbr.append(np.exp(-(dist[indnbr[p]]**2)/np.sqrt(epsilon)))
    return(dnbr,wnbr,indnbr,indtri)

def findbase(Da,dnbr,wnbr,indnbr):
    Oplist = []
    for p in range(len(Da)):
        if len(indnbr[p]) < 2:
            Oplist.append(np.array([-1]))
            continue
        Xp = (Da[indnbr[p]] - Da[p]).T # Only near neighbor is used 
        Dp = np.diag(np.sqrt(wnbr[p]))    
        Bp = np.dot(Xp, Dp)
        u, sigma, v = npl.svd(Bp)
        Op = u[:,:dim]
        Oplist.append(Op)
    return(Oplist)

def angSum(i,j,k,Oplist,Da):
    basei = Oplist[i].T
    basej = Oplist[j].T
    basek = Oplist[k].T
    if((basei==-1).all() or (basej == -1).all() or (basek == -1).all()):
        return(-1)
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
    return(angi+angj+angk,iij,iik)
    
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
 
def curvatp_gb_simp(p,ind,Oplist,dnbr):
    if(p in ind):
        ind = np.setdiff1d(ind,np.array(p))
        ang = angSum(p,ind[0],ind[1],Oplist,Da)
        Area = findArea(p,ind[0],ind[1],dnbr)
        if(Area < 0 or ang[0] < 0):
            return(None)
        else: 
            res = (ang[0]-np.pi)/Area
    else:
        ang = angSum(ind[0],ind[1],ind[2],Oplist,Da)
        Area = findArea(p,ind[0],ind[1],dnbr)+findArea(p,ind[1],ind[2],dnbr)+findArea(p,ind[0],ind[2],dnbr)
        if(Area < 0 or ang[0] < 0):
            print("Exit Code 2")
            return(None)
        else: 
            res = (ang[0]-np.pi)/Area
    coef = np.kron(ang[1],np.kron(ang[2],np.kron(ang[2],ang[1])))/(npl.norm(ang[1])**2*npl.norm(ang[2])**2-(np.inner(ang[1],ang[2]))**2)
    return(res,coef)
    
def psihuber(u,M=1.35):
    if(abs(u) <= M):
        return(u)
    elif(u > M):
        return(M)
    else:
        return(-M)
    
def huberfit(y,x,maxiter=1000,tol=1e-8):
    iteration = 1
    y=y.T
    try:
        beta = npl.solve(np.dot(x.T,x)+0.001*np.eye(x.shape[1]),np.dot(x.T,y)).T
    except npl.linalg.LinAlgError as err:
        print("Exit Code 3")
        return(None)
    for i in range(maxiter):
        beta_prev = beta
        resi = y - np.dot(x,beta)
        mad = np.median(abs(resi - np.median(resi)))
        scale = mad/0.6745
        weight = np.diag(np.array(map(psihuber,resi/scale))/(resi/scale))
        temp = np.dot(x.T,weight)
        try:
            beta = npl.solve(np.dot(temp,x)+0.001*np.eye(x.shape[1]),np.dot(temp,y))
        except npl.linalg.LinAlgError as err:
            print("Exit Code 4")
            return(None)
        if(sum(abs(beta_prev-beta))<tol):
            break
    return(beta)
    
def curvatp_gb_ensemble(p,indtri,Da,Oplist,dnbr):
    nbr = indtri[p]
    length = len(nbr)
    if(length < 2):
        print("Exit Code 5")
        return(None)
    sample_res = []
    sample_coef = []
    for indi in range(length):
        for indj in range(indi+1,length):
            ind = np.array([p,nbr[indi],nbr[indj]])
            temp = curvatp_gb_simp(p,ind,Oplist,dnbr)
            if(temp is None):
                continue
            sample_res.append(temp[0])
            sample_coef.append(temp[1])
    sample_res=np.array(sample_res)
    sample_coef = np.array(sample_coef)
    if(np.sum(~np.isnan(sample_res))==0):
        print("Exit Code 6")
        return(None)
    else:
        W = np.dot(sample_coef,S)
#       Finding the Huber result
        res = huberfit(sample_res,W)
        print(res)
    return(res)

t1 = time.time()  
res_huber = []
for ind in range(repeat_time):
    Da = gendata(mean,cov,size)
    num_of_nbr = np.zeros(len(Da))
    dnbr,wnbr,indnbr,indtri = findnbr(Da,epsilon)
    Oplist = findbase(Da,dnbr,wnbr,indnbr)
    temp = curvatp_gb_ensemble(0,indtri,Da,Oplist,dnbr)
    if(temp is not None):
        res_huber.append(temp)
res_huber = np.vstack(res_huber)
for i in range(6):
    print(np.median(res_huber[:,i]))
t2 = time.time()
print("time: ",t2-t1)


    

