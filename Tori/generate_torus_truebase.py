#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
st = '''#!/bin/bash
#SBATCH --job-name torus_{0}_{1}
#SBATCH --partition short
#SBATCH --ntasks 1
#SBATCH --time 0-11:00
#SBATCH --mem-per-cpu=30000
#SBATCH -o output_%j.out
#SBATCH -e error_%j.err
module load python
echo “Running plot script on a single CPU core”
python torus_{0}_{1}.py
'''

with open('torus.py') as f:
    fr = f.read()
    
size = [5000,10000,15000]
eps = [0.3,0.4,0.5,0.6]
for i in size:
    for j in eps:
        with open('torus_'+str(i)+'_'+str(int(10*j))+'.py','w') as f:
            f.write(fr.format(i,j))
        with open('torus_'+str(i)+'_'+str(int(10*j))+'.sbatch','w') as fi:
            fi.write(st.format(i,int(10*j)))
        os.spawnlp(os.P_NOWAIT,'sbatch','sbatch','torus_'+str(i)+'_'+str(int(10*j))+'.sbatch')
