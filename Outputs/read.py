# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 09:56:44 2017

@author: hanyuz6
"""

import os
import re
import numpy as np
my_path = os.getcwd()
'''
dat_sphere4d = []
for id in range(177,253,1)+range(255,277,1):
    print(id)
    with open(my_path+'/output_60'+str(id)+'.out','r') as f:
        cont = f.readlines()
    pat = re.compile(r"([-+]?\d*\.\d*(e[-+]?\d*)?)")
    dat_temp = []
    for i in range(1,len(cont)-19):
        temp = re.findall(pat,cont[i])
        for k in range(len(temp)):
            res = float(temp[k][0])    
            dat_temp.append(res)
    dat_temp = np.reshape(np.asarray(dat_temp),(100,20))
    dat_sphere4d.append(dat_temp)
dat_sphere4d = np.vstack(dat_sphere4d)

dat_torus = []
for id in range(61066,61116,1):
    print(id)
    with open(my_path+'/output_'+str(id)+'.out','r') as f:
        cont = f.readlines()
    pat = re.compile(r"([-+]?\d*\.\d*(e[-+]?\d*)?)")
    dat_temp = []
    for i in range(1,len(cont)-10):
        temp = re.findall(pat,cont[i])
        for k in range(len(temp)):
            res = float(temp[k][0])    
            dat_temp.append(res)
    dat_temp = np.reshape(np.asarray(dat_temp),(100,6))
    dat_torus.append(dat_temp)
dat_torus = np.vstack(dat_torus)

dat_ctorus = []
for id in range(62608,62707,1):
    print(id)
    with open(my_path+'/output_'+str(id)+'.out','r') as f:
        cont = f.readlines()
    pat = re.compile(r"([-+]?\d*\.\d*(e[-+]?\d*)?)")
    dat_temp = []
    for i in range(len(cont)-25,len(cont)-1):
        temp = re.findall(pat,cont[i])
        for k in range(len(temp)):
            res = float(temp[k][0])    
            dat_temp.append(res)
    dat_temp = np.asarray(dat_temp)
    dat_ctorus.append(dat_temp)
dat_ctorus = np.hstack(dat_ctorus)
'''
dat_psphere = []
for id in range(64793,64874,1)+range(64875,64893,1):
    print(id)
    with open(my_path+'/output_'+str(id)+'.out','r') as f:
        cont = f.readlines()
    pat = re.compile(r"([-+]?\d*\.\d*(e[-+]?\d*)?)")
    dat_temp = []
    for i in range(len(cont)-25,len(cont)-1):
        temp = re.findall(pat,cont[i])
        for k in range(len(temp)):
            res = float(temp[k][0])    
            dat_temp.append(res)
    dat_temp = np.asarray(dat_temp)
    dat_psphere.append(dat_psphere)
dat_psphere = np.hstack(dat_psphere)
