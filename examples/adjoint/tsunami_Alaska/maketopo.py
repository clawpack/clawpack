"""
Download topo and dtopo files needed for this example.
    
Call functions with makeplots==True to create plots of topo, slip, and dtopo.
"""

import os

import clawpack.clawutil.data

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')

def get_topo(makeplots=False):
    """
    Retrieve the topo file from online.
    """
    from clawpack.geoclaw import topotools
    topo_fname = 'etopo1min170E124W40N61N.asc'
    url = 'http://students.washington.edu/bndavis/misc/topo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir, 
            file_name=topo_fname, verbose=True)
            
    topo_fname = 'etopo4min120E110W0N62N.asc'
    url = 'http://students.washington.edu/bndavis/misc/topo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir,
            file_name=topo_fname, verbose=True)
            
    topo_fname = 'cc-1sec-c.asc'
    url = 'http://students.washington.edu/bndavis/misc/topo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir,
            file_name=topo_fname, verbose=True)
            
    topo_fname = 'cc-1_3sec-c_pierless.asc'
    url = 'http://students.washington.edu/bndavis/misc/topo/' + topo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir,
            file_name=topo_fname, verbose=True)

    if makeplots:
        from matplotlib import pyplot as plt
        topo = topotools.Topography(topo_fname, topo_type=2)
        topo.plot()
        fname = os.path.splitext(topo_fname)[0] + '.png'
        plt.savefig(fname)
        print "Created ",fname
    
def make_dtopo(makeplots=False):
    """
    Create dtopo data file for deformation of sea floor due to earthquake.
    Uses the Okada model with fault parameters and mesh specified below.
    """
    from clawpack.geoclaw import dtopotools
    import numpy

    dtopo_fname = 'AASZ04v2.tt3'
    url = 'http://students.washington.edu/bndavis/misc/dtopo/alaska/' + dtopo_fname
    clawpack.clawutil.data.get_remote_file(url, output_dir=scratch_dir,
                                           file_name=dtopo_fname, verbose=True)

    if makeplots:
        from matplotlib import pyplot as plt
        if fault.dtopo is None:
            # read in the pre-existing file:
            print "Reading in dtopo file..."
            dtopo = dtopotools.DTopography()
            dtopo.read(dtopo_fname, dtopo_type=3)
            x = dtopo.x
            y = dtopo.y
        plt.figure(figsize=(12,7))
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(122)
        fault.plot_subfaults(axes=ax1,slip_color=True)
        ax1.set_xlim(x.min(),x.max())
        ax1.set_ylim(y.min(),y.max())
        dtopo.plot_dz_colors(1.,axes=ax2)
        fname = os.path.splitext(dtopo_fname)[0] + '.png'
        plt.savefig(fname)
        print "Created ",fname


if __name__=='__main__':
    get_topo(False)
    make_dtopo(False)
