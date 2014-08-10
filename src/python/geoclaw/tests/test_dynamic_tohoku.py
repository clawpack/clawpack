
from clawpack.geoclaw import dtopotools
from numpy import linspace, array, load
import matplotlib.pyplot as plt
import time
import os

datadir = os.path.join(os.environ['CLAW'], 'geoclaw', 'datafiles')

shoreline_fname = os.path.join(datadir,'tohoku_shoreline_1min.npy')
shoreline_xy = load(shoreline_fname)

subfault_fname = os.path.join(datadir,'tohoku_ucsb.txt')
fault = dtopotools.UCSBFault()
fault.read(subfault_fname)
fault.rupture_type = 'dynamic'

fault.plot_subfaults(slip_color=True)  # plot final slip
plt.show()

# seafloor deformation:

xlower = 140.
xupper = 146.
ylower = 35.
yupper = 41.
xylim = [xlower,xupper,ylower,yupper]

# dtopo parameters for 4 min resolution:
mx = int((xupper - xlower)*60 + 1)
my = int((yupper - ylower)*60 + 1)


x = linspace(xlower,xupper,mx)
y = linspace(ylower,yupper,my)

tmax = 0.
for s in fault.subfaults:
    tmax = max(tmax, s.rupture_time + s.rise_time + s.rise_time_ending)
print "rupture ends at time ",tmax

times = linspace(0,tmax,10)
dtopo = fault.create_dtopography(x,y,times,verbose=True)

dz_final = dtopo.dz_list[-1]
dz_max = dz_final.max()

dtopo.plot_dz_colors(t=tmax)
plt.show()

# Incorporate this function in dtopotools to replace animate_dz_colors?
def plot_subfaults_dz(t, fig=None):
    if fig is None:
        fig = plt.figure(figsize=(12,5))
    else:
        fig.clf()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    fault.plot_subfaults(axes=ax1, slip_color=True, slip_time=t, xylim=xylim)
    dtopo.plot_dz_colors(axes=ax2, t=t, cmax_dz=dz_max)
    ax1.plot(shoreline_xy[:,0],shoreline_xy[:,1],'g')
    ax2.plot(shoreline_xy[:,0],shoreline_xy[:,1],'g')
    plt.axis(xylim)
    fig.show()
    return fig
    
fig = plt.figure(figsize=(12,5))

for t in list(linspace(0,150,16)) + [170,200]:
    plot_subfaults_dz(t,fig)
    plt.draw()
    plt.show()
    time.sleep(1)

fname = 'tohoku_ucsb_dynamic.tt3'
dtopo.write(fname, 3)
print 'Created ',fname


