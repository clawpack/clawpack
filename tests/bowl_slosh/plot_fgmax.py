
"""
Plot fgmax output from GeoClaw run.

"""

from __future__ import absolute_import
from pylab import *
from numpy import ma
from clawpack.geoclaw import fgmax_tools

fg = fgmax_tools.FGmaxGrid()
fg.read_input_data('fgmax1.txt')
fg.read_output()

figure(1)
clf()
surface = ma.masked_where(fg.h < 0.001, fg.h + fg.B)
contourf(fg.X,fg.Y,surface,10)
cb = colorbar()
cb.set_label('meters')
title('Max surface elevation')

figure(2)
clf()
#s = ma.masked_where(fg.s<-1e10, fg.s)
contourf(fg.X,fg.Y,fg.s,10)
cb = colorbar()
cb.set_label('meters / sec')
title('Max speed')

figure(3)
clf()
contourf(fg.X,fg.Y,fg.arrival_time,10)
cb = colorbar()
cb.set_label('seconds')
title('Arrival time')

