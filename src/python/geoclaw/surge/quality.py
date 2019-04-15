"""Module for QA/QC on GEOCLAW model outputs."""

from glob import glob
from os.path import join
from clawpack.pyclaw import Solution
import numpy as np


class InstabilityError(Exception):
    pass

def check_model_stability(model_output_dir,
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

    return None