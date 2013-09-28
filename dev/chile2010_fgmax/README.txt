
Test case for fgmax routines.  Needs to be cleaned up and moved to apps
repository.

To test:

python make_fgmax.py
make .output
python plot_fgmax_grid.py
python plot_fgmax_transect.py

This should produce the files 
   zeta.png             maximum amplitude along with contours of arrival times
   arrival_times.png    color map of arrival times
   zetatimes.png        color map of time of maximum amplitude
   zeta_transect.png    1d plot of solution on a transect (from FG2)
   
