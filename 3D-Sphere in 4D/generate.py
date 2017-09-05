#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import numpy as np
st = '''#!/bin/bash
#SBATCH —job-name 4dsphere_{0}_{1}
#SBATCH —partition short
#SBATCH —ntasks 1
#SBATCH —time 0-11:00
#SBATCH —mem-per-cpu=30000
#SBATCH -o output_%j.out
#SBATCH -e error_%j.err
module load python
echo “Running plot script on a single CPU core”
python sphere4d_{0}_{1}.py
'''

np.set_printoptions(precision = 6, threshold=np.inf)
with open('sphere4d.py') as f:
    fr = f.read()
    
size = [1000,5000,10000]
eps = [0.1,0.2,0.3,0.4,0.5]
for i in size:
    for j in eps:
        with open('sphere4d_'+str(i)+'_'+str(int(10*j))+'.py','w') as f:
            f.write(fr.format(i,j))
        with open('sphere4d_'+str(i)+'_'+str(int(10*j))+'.sbatch','w') as fi:
            fi.write(st.format(i,int(10*j))
        os.spawnlp(os.P_NOWAIT,'sbatch','sbatch','sphere4d_'+str(i)+'_'+str(int(10*j))+'.sbatch')