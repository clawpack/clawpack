#!/usr/bin/env python
# encoding: utf-8

import numpy as np

plot = False

raw_bathy_file = 'Survey_Data_SWG/SWG_surveys.xyz'

data = np.loadtxt(raw_bathy_file,skiprows=1)

x = data[:,0]
y = data[:,1]
z = data[:,2]

xll = x[0]
yll = y[0]
num_cells = [0,0]
for coord in x[1:]:
    if xll == coord:
        print coord,xll
        break
    num_cells[0] += 1
    
print num_cells
        

