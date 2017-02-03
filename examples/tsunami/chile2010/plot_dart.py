"""
Compare DART buoy results to NOAA results.
"""

from __future__ import absolute_import
import pylab 
from pyclaw.plotters.data import ClawPlotData
from matplotlib import image

dartpng = image.imread('dart32412_comp-2.png')
pylab.figure(400)
pylab.clf()
pylab.imshow(dartpng)
pylab.hold(True)

# origin of plot in pixels on image:
t0 = 195.
y0 = 221.
tscale = 120.0
yscale = -446.66666666666669
plotdata = ClawPlotData()
plotdata.outdir = "_output"
gaugedata = plotdata.getgauge(32412)
t = gaugedata.t
eta = gaugedata.q[:,3]
t = t0 + tscale * (t / 3600.)
y = y0 + yscale*eta

pylab.plot(t,y,'b')
pylab.plot([160,225],[550,550],'b')
pylab.text(240,560,'GeoClaw (added to NOAA original)',fontsize=8)
pylab.text(100,600,'NOAA Original from http://nctr.pmel.noaa.gov/chile20100227/dart32412_comp-2.pdf',fontsize=8)


pylab.xlim([50,1300])
pylab.ylim([640,0])
pylab.axis('off')

pylab.savefig('dart.png')
