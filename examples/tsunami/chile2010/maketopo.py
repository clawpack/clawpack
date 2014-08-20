"""
Create topo files needed for this example:
    etopo10min120W60W60S0S.asc        download from GeoClaw topo repository
    usgs100227.tt1                    create using Okada model and .cfg file
    
"""

import os,sys

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'SCRATCH')

def get_topo():
    """
    Retrieve the topo file from the GeoClaw repository.
    """
    from clawpack.geoclaw import topotools, util
    topo_fname = 'etopo10min120W60W60S0S.asc'
    url = 'http://www.geoclaw.org/topo/etopo/' + topo_fname
    util.get_remote_file(url, output_dir=scratch_dir, file_name=topo_fname)

    
def make_dtopo():
    """
    Create dtopo data file for deformation of sea floor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified below.
    """
    from clawpack.geoclaw import dtopotools
    import numpy

    dtopo_fname = os.path.join(scratch_dir, "dtopo_usgs100227.tt3")

    if os.path.exists(dtopo_fname):
        print "*** Not regenerating dtopo file (already exists): %s" % dtopo_fname
    else:
        print "Using Okada model to create dtopo file"

        usgs_subfault = dtopotools.SubFault()
        usgs_subfault.strike = 16.
        usgs_subfault.length = 450.e3
        usgs_subfault.width = 100.e3
        usgs_subfault.depth = 35.e3
        usgs_subfault.slip = 15.
        usgs_subfault.rake = 104.
        usgs_subfault.dip = 14.
        usgs_subfault.longitude = -72.668
        usgs_subfault.latitude = -35.826
        usgs_subfault.coordinate_specification = "top center"

        x = numpy.linspace(-77, -67, 100)
        y = numpy.linspace(-40, -30, 100)
        times = [0., 1.]

        fault = dtopotools.Fault()
        fault.subfaults = [usgs_subfault]
        fault.create_dtopography(x,y,times)
        dtopo = fault.dtopo
        dtopo.write(dtopo_fname, 3)

        print "Mw = ",fault.Mw()

        if 1:
            from matplotlib import pyplot as plt
            plt.figure(figsize=(12,7))
            ax1 = plt.subplot(121)
            ax2 = plt.subplot(122)
            fault.plot_subfaults(axes=ax1,slip_color=True)
            ax1.set_xlim(x.min(),x.max())
            ax1.set_ylim(y.min(),y.max())
            dtopo.plot_dz_colors(1.,axes=ax2)
            fname = 'dtopo_usgs100227.png'
            plt.savefig(fname)
            print "Created ",fname


if __name__=='__main__':
    get_topo()
    make_dtopo()
