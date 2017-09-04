#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 09:58:12 2017

@author: mac
"""

import os
import numpy as np
st = '''
#!/bin/bash
#SBATCH —job-name size{0}_{1}
#SBATCH —partition short
#SBATCH —ntasks 1
#SBATCH —time 0-10:00
#SBATCH —mem-per-cpu=30000
#SBATCH -o output_%j.out
#SBATCH -e error_%j.err
module load python
echo “Running plot script on a single CPU core”
python Sim_gb2_{0}_{1}.py
'''


with open('Sim_gb2.py') as f:
    fr = f.read()
    
size = [1000,5000,10000,15000]
eps = 0.1 + np.arange(11) * 0.025
for i in size:
    for j in eps:
        with open('Sim_gb_'+str(i)+'_'+str(1000*j)+'.py','w') as f:
            f.write(fr.format(i,j))
        with open('Sim_gb_'+str(i)+'_'+str(1000*j)+'.sbatch','w') as fi:
            fi.write(st.format(i,1000*j))
 #       os.spawnlp(os.P_NOWAIT,'sbatch','sbatch','Sim_gb_'+str(i)+'_'+str(1000*j)+'.sbatch')