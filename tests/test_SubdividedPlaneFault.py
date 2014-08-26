
from clawpack.geoclaw import dtopotools
from matplotlib import pyplot as plt

# get a unit source fault plane as starting point:
sift_slip = {'acsza1':1.}
fault = dtopotools.SiftFault(sift_slip)
fault_plane = fault.subfaults[0]
Mo = fault_plane.Mo()
print "original Mo = ",Mo

fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3)
print "new Mo = ",fault2.Mo()
fault2.plot_subfaults(slip_color=True)
plt.show()

slip_function = lambda xi,eta: xi*(1-xi)*eta
fault2 = dtopotools.SubdividedPlaneFault(fault_plane, nstrike=5, ndip=3, 
        slip_function=slip_function, Mo=Mo)
print "new Mo = ",fault2.Mo()
fault2.plot_subfaults(slip_color=True)
plt.show()

fault2.subdivide(nstrike=20, ndip = 10, slip_function=slip_function, Mo=Mo)
print "with finer resolution, Mo = ",fault2.Mo()
fault2.plot_subfaults(slip_color=True)
plt.show()
