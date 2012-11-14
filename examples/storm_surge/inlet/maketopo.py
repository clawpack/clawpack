#!/usr/bin/env python
# encoding: utf-8
""" 
Functions for creating topography
    
"""

import os
import sys
import numpy as np

from clawpack.clawutil.data import Data
import geoclaw.topotools as tt

# ============================================================================
#  Profiles
def generate_profiles():
    # Parameters
    data = Data('./amr2ez.data')
    dx = abs(data.xupper-data.xlower) / (data.mx)
    dy = abs(data.yupper-data.ylower) / (data.my)
    d = min(dx,dy)
    mx = int((data.xupper-data.xlower) / d) + 8
    my = int((data.yupper-data.ylower) / d) + 8

    xlower = data.xlower-d*4.0
    ylower = data.ylower-d*4.0
    xupper = data.xupper+d*4.0
    yupper = data.yupper+d*4.0

    beach_slope = 0.05
    basin_depth = -3000
    shelf_depth = -200
    x0 = 350e3
    x1 = 450e3
    x2 = 480e3

    shelf_slope = (basin_depth - shelf_depth) / (x0 - x1)
    y_end = beach_slope * (xupper - 477e3) - 100.0
    
    # Stommel problem depth
    stommel_depth = -1000
    
    # Points and depths for gulf shelf
    #  1: (25°39'2.85"N, 86° 7'24.77"W)   --   -3228 m   --   0.0 m
    #  2: (27°53'44.74"N, 88° 0'34.02"W)   --   -2438 m   --   312.17313 km
    #  3: (28°59'47.14"N, 88°59'53.19"W)   --   -188 m   --    467.59957 km
    #  4: ( 29° 4'6.90"N,  89° 4'11.39"W)    --   0 m   --   479.10557 km

    bathy_profiles = {"simple_shelf":[(xlower,basin_depth),
                                      (x0,basin_depth),
                                      (x1,shelf_depth),
                                      (x2,shelf_depth),
                                      (xupper,beach_slope*(xupper-x2)+shelf_depth)],
                      "shallow_shelf":[(xlower,-100),
                                       (x2,-100),
                                       (xupper,y_end)],
                      "gulf_shelf":[(0.0,-3228),
                                    (312e3,-2438),
                                    (467e3,-188),
                                    (479e3,0.0),
                                    (579e3,300.0)],
                      "step_shelf1":[(0.0,-2000.0),
                                     (470e3-0.001,-2000.0),
                                     (470e3,-200.0),
                                     (500e3,-200.0)],
                      "shelf":[(0.0,-3000),
                               (400e3,-2700),
                               (450e3,-100),
                               (500e3,-100)],
                      "continental_shelf":[(2000e3,-7000),
                                           (2800e3,-3000),
                                           (2900e3,-100),
                                           (3000e3,0.0)],
                      "flat100":[(0.0,-100),
                              (1000e3,-100)],
                      "flat1000":[(0.0,-1000),
                              (1000e3,-1000)],
                      "flat_stommel":[(0.0,stommel_depth),
                                      (1000e3,stommel_depth)]}
    return bathy_profiles

# ============================================================================
#  Topography generation functions 
# ============================================================================
def write_topo_file(topo_file,topo_type=1,factor=4,bathy_type=None,plot=False,
                        force=False,verbose=False):
    """Creates topography file needed by the simulation"""
    
    if os.path.exists(topo_file):
        if not force:
            ans = raw_input("The file %s already exists, replace: [Y/n] " % topo_file)
            if ans[0].capitalize() == 'N':
                print "File %s already exists, use force if you want to replace the existing file." % topo_file
                return
    print "Creating topography file ",topo_file
    
    # Look at parameters of domain
    data = Data('./amr2ez.data')
    dx = abs(data.xupper-data.xlower) / (data.mx)
    dy = abs(data.yupper-data.ylower) / (data.my)
    d = min(dx,dy)
    mx = int((data.xupper-data.xlower) / d) + 8
    my = int((data.yupper-data.ylower) / d) + 8

    xlower = data.xlower-d*4.0
    ylower = data.ylower-d*4.0
    xupper = data.xupper+d*4.0
    yupper = data.yupper+d*4.0
    
    bathy_profiles = generate_profiles()

    # Pick out bathymetry profile
    if bathy_type is None:
        bathy_profile = 'flat'
    if bathy_type in bathy_profiles.keys():
        bathy_profile = bathy_profiles[bathy_type]
    else:
        err_msg = "Bathy profile %s not found, available options are" % bathy_type
        for name in bathy_profiles.keys():
            err_msg = '\n   '.join((err_msg,name))
        raise Exception(err_msg)
    if verbose:
        print "%s profile: %s" % (bathy_type,bathy_profile)

    # Create and write topography
    bathy_func = tt.create_topo_func(bathy_profile,verbose=verbose)
    N = len(bathy_profile)
    if topo_type == 1:
        tt.topo1writer(topo_file,bathy_func,xlower,xupper,ylower,yupper,
            factor*mx,factor*my)
    elif topo_type == 2:
        print (xlower - xupper) / (factor*mx)
        print (ylower - yupper) / (factor*my)
        tt.topo2writer(topo_file,bathy_func,xlower,xupper,ylower,yupper,
            factor*mx,factor*my)
    else:
        raise Exception("Unknown topo type %s requested" % topo_type)

    # Plotting
    if plot:
        plot_profiles({bathy_type:bathy_profile},factor=factor)
        
        
def plot_profiles(profiles,topo_type=1,factor=4,verbose=True):
    
    import tempfile
    import shutil
    import matplotlib.pyplot as plt
    
    # Assume an amr2ez.data file is present
    data = Data('./amr2ez.data')
    dx = abs(data.xupper-data.xlower) / (data.mx)
    dy = abs(data.yupper-data.ylower) / (data.my)
    d = min(dx,dy)
    mx = int((data.xupper-data.xlower) / d) + 8
    my = int((data.yupper-data.ylower) / d) + 8

    xlower = data.xlower-d*4.0
    ylower = data.ylower-d*4.0
    xupper = data.xupper+d*4.0
    yupper = data.yupper+d*4.0
    
    bathy_path = tempfile.mkdtemp()
    grid_size = int(np.ceil(np.sqrt(len(profiles.keys()))))
    index = 1
    plt.figure(1)
    
    for (name,profile) in profiles.iteritems():
        topo_file = os.path.join(bathy_path,name)
        
        # Write out temporary file
        bathy_func = tt.create_topo_func(profile,verbose=verbose)
        N = len(profile)
        factor = 4
        if topo_type == 1:
            tt.topo1writer(topo_file,bathy_func,xlower,xupper,ylower,yupper,
                        factor*mx,factor*my)
        elif topo_type == 2:
            tt.topo2writer(topo_file,bathy_func,xlower,xupper,ylower,yupper,
                        factor*mx,factor*my)
        else:
            raise Exception("Invalid topo_type %s." % topo_type)
                        
        # Read it back in
        (X,Y,Z) = tt.topofile2griddata(topo_file,topotype=topo_type)
        slice_index = int(X.shape[0] / 2.0)
        x = X[slice_index,:]
        b = Z[slice_index,:]
        
        limits = [np.min(x),np.max(x),np.min(b)*1.05,np.max(b)*1.05]
        wet_index = np.nonzero(b < 0.0)
        
        plt.subplot(grid_size,grid_size,index)
        plt.hold(True)
        plt.plot(x,b,'k',linewidth=2)
        plt.plot(x[wet_index],np.zeros(x[wet_index].shape),'b')
        plt.axis(limits)
        plt.title(name)
        plt.xlabel('km')
        plt.ylabel('m')
        locs,labels = plt.xticks()
        labels = locs/1e3
        plt.xticks(locs,labels)
        plt.hold(False)
        index = index + 1
        
    plt.show()
        
    shutil.rmtree(bathy_path)

        
if __name__ == "__main__":
    # plot_profiles(generate_profiles())
    # write_bathy_profile('./topo.data',None,plot=True,force=True)
    # Create bathy
    topo_file = './topo.tt1'
    if len(sys.argv) > 1:
        topo_file = sys.argv[1]
    write_topo_file(topo_file,topo_type=1,factor=4,
                            bathy_type='simple_shelf',plot=False,force=True)
