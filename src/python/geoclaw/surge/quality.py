"""Module for QA/QC on GEOCLAW model outputs."""

from glob import glob
from os.path import join
from clawpack.pyclaw import Solution
import numpy as np


class InstabilityError(Exception):
    pass

def quality_check(model_output_dir,
                  regions_to_check = [[-81,22,-40,55], # atlantic
                                      [-100,15,-82,32], # gulf
                                      [-88.5,13.25,-70,19.75]], # carribean
                  frames_to_check = 4,
                  mean_tol_abs = .01,
                  min_depth = 300):
    """Run a simple check to flag runs that were potentially numerically unstable.
    This function looks through the final *frames_to_check* frames of a model output
    and checks all cells in AMR level 1 that are below *min_depth* and 
    fall within each region in the *regions_to_check* list (which are intended 
    to encompass individual basins). If the absolute value of the mean surface height
    for these cells, which should not really have much surge, is above/mean_tol_abs, 
    it raises an error.
    
    :Input:
    - *model_output_dir* (str) Path to the output directory for a given model
    - *regions_to_check* (list of lists of floats) A list of 4-element lists
        which define the (left, bottom, right, top) of regions to check within.
    - *frames_to_check* (int) The number of final frames of the model output to check
    - *mean_tol_abs* (float) Model is flagged if the absolute value of the mean surface
        height within any region in *regions_to_check* is above this value (meters).
    - *min_depth* (float) Only cells that are below *min_depth* (meters) are checked, in 
        order to avoid including cells that potentially should have large surge values.
        
    :Raises:
    - InstabilityError if the model is flagged as potentially unstable in one of the regions.
    """
    
    # find which frame is last
    output_files = glob(join(model_output_dir,'fort.b*'))
    last_frame = int(output_files[-1].split('/')[-1][-4:])

    # get sl_init
    with open(join(model_output_dir,'geoclaw.data'),'r') as f:
        for l in f:
            if '=: sea_level' in l:
                sl_init = float(l.strip().split()[0])

    # bring in solution
    soln = Solution()
    for frame in range(last_frame-frames_to_check+1, last_frame+1):
        soln.read(frame, path = model_output_dir, file_format='binary', read_aux=False)

        for r in regions_to_check:
            all_vals = np.array([])
            for s in soln.states:

                # only looking at lowest AMR level
                if s.patch.level >1:
                    continue

                x = s.grid.c_centers[0]
                y = s.grid.c_centers[1]
                eta = s.q[3,:,:]
                topo = eta - s.q[0,:,:]
                
                # only count 
                eta = np.where(topo<(-min_depth),eta,np.nan)
                in_reg = eta[np.ix_((x[:,0]>=r[0]) & (x[:,0]<=r[2]),
                             (y[0,:]>=r[1]) & (y[0,:]<=r[3]))]
                all_vals = np.concatenate([all_vals,in_reg.flatten()])

            # adjust for sl_init
            all_vals = all_vals - sl_init

            if all_vals.shape[0]>0:
                if abs(np.nanmean(all_vals)) > mean_tol_abs:
                    raise InstabilityError("Model possibly unstable due to large magnitude deep "
                                           "ocean surge at end of run ({:.1f} cm in region {})".format(
                                           np.nanmean(all_vals)*100, r))

                    
def get_max_boundary_fluxes(model_output_dir, max_refinement_depth_to_check = None, 
                            frames_to_check = 1):
    """A common cause of instability is a persistant normal flux at the boundary.
    This code checks the last frame(s) of a model output and returns the maximum
    magnitude normal fluxes and currents at each boundary, as well as the maximum 
    magnitude fluxes and currents observed within the entire domain.
    
    :Input:
    - *model_output_dir* (str) Path to the output directory for a given model
    - *max_refinement_depth_to_check* (int or None, optional) How many refinement levels to 
        loop through to find max values. Runs quicker when just checking level 1, but
        you may find higher max values at higher refinement levels. None (default)
        means check all levels
    - *frames_to_check* (int) How many of the last output frames to check for max. Default
        is just one.
    """

    # find which frame is last
    output_files = glob(join(model_output_dir,'fort.q*'))
    frames = [int(i.split('/')[-1][-4:]) for i in output_files]
    last_frame = max(frames)
    
    # get domain and file format
    with open(join(model_output_dir,'claw.data'),'r') as f:
        for l in f:
            if '=: lower' in l:
                xl,yl = [float(i) for i in l.strip().split()[:2]]
            elif '=: upper' in l:
                xu,yu = [float(i) for i in l.strip().split()[:2]]
            elif '=: output_format' in l:
                of = int(l.strip().split()[0]) - 1
    file_format = ['ascii','netcdf','binary'][of]

    maxhxl = -np.inf
    maxhyl = -np.inf
    maxcurrxl = -np.inf
    maxcurryl = -np.inf
    maxhxu = -np.inf
    maxhyu = -np.inf
    maxcurrxu = -np.inf
    maxcurryu = -np.inf
    maxhx_overall = -np.inf
    maxhy_overall = -np.inf
    maxcurrx_overall = -np.inf
    maxcurry_overall = -np.inf

    for f in range(last_frame+1-frames_to_check,last_frame+1):

        soln = Solution()
        soln.read(f, path = model_output_dir, file_format=file_format, read_aux=False)

        for s in soln.states:

            # only looking at lowest AMR levels
            if max_refinement_depth_to_check is not None:
                if s.patch.level > max_refinement_depth_to_check:
                    continue
                    
            # get rounding error tolerance
            delta = s.grid.dimensions[0].delta
            edge_tol = delta * .001
            
            x = s.grid.c_centers[0]
            xedges = s.grid.c_nodes[0]
            y = s.grid.c_centers[1]
            yedges = s.grid.c_nodes[1]

            eta = s.q[3,:,:]
            h = s.q[0,:,:]
            hx = s.q[1,:,:]
            curr_x = hx / h
            hy = s.q[2,:,:]
            curr_y = hy / h
            topo = eta - s.q[0,:,:]
            maxhx_overall = np.nanmax([np.abs(hx).max(),maxhx_overall])
            maxhy_overall = np.nanmax([np.abs(hy).max(),maxhy_overall])
            maxcurrx_overall = np.nanmax([np.nanmax(np.abs(curr_x)),maxcurrx_overall])
            maxcurry_overall = np.nanmax([np.nanmax(np.abs(curr_y)),maxcurry_overall])

            if abs(xedges[0,0] - xl) < edge_tol:
                maxhxl = np.nanmax([maxhxl,np.nanmax(np.abs(hx[0,:]))])
                maxcurrxl = np.nanmax([maxcurrxl,np.nanmax(np.abs(curr_x[0,:]))])
            if abs(xedges[-1,0] - xu) < edge_tol:
                maxhxu = np.nanmax([maxhxu,np.nanmax(np.abs(hx[-1,:]))])
                maxcurrxu = np.nanmax([maxcurrxu,np.nanmax(np.abs(curr_x[-1,:]))])
            if abs(yedges[0,0] - yl) < edge_tol:
                maxhyl = np.nanmax([maxhyl,np.nanmax(np.abs(hy[:,0]))])
                maxcurryl = np.nanmax([maxcurryl,np.nanmax(np.abs(curr_y[:,0]))])
            if abs(yedges[0,-1] - yu) < edge_tol:
                maxhyu = np.nanmax([maxhyu,np.nanmax(np.abs(hy[:,-1]))])
                maxcurryu = np.nanmax([maxcurryu,np.nanmax(np.abs(curr_y[:,-1]))])
                
    return ({'max_normal_fluxes':{'W': maxhxl,
                              'E': maxhxu,
                              'N': maxhyu,
                              'S': maxhyl},
             'max_normal_currents':{'W': maxcurrxl,
                                'E': maxcurrxu,
                                'N': maxcurryu,
                                'S': maxcurryl},
             'domain_maxs': {'hu': maxhx_overall,
                             'hv': maxhy_overall,
                             'u': maxcurrx_overall,
                             'v': maxcurry_overall}})