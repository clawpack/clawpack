
import maketopo
from pylab import *

x=linspace(0.98,1.04,121) 
y=linspace(0.96,1.03,141)
X,Y = meshgrid(x,y)
Z = maketopo.topo(X,Y)
contour(X,Y,Z,linspace(-10,30,41),colors='k')
plot([1.008],[1],'bo')
plot([1.01],[1],'ro')
plot([1.012],[1],'go')
title('1 meter contours and 3 gauges')
savefig('bay.png')
print "Created bay.png"

