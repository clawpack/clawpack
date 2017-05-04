
from __future__ import absolute_import
from pylab import *

outdir = '_output'
d = loadtxt(outdir + '/fort.FG1.valuemax')
b = loadtxt(outdir + '/fort.FG1.aux1')
topo = b[:,4]  # for level 3
figure(401)
clf()
