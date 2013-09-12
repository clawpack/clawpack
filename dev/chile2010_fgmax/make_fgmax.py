"""
Create fgmax_grid.txt input file 
"""

from numpy import linspace


def make_fgmax_grid():
    xlower = -120
    xupper = -60
    ylower = -60
    yupper = 0
    dx = 0.2
    dy = dx
    mx = int(round((xupper-xlower)/dx)) + 1  
    my = int(round((yupper-ylower)/dy)) + 1  
    print "Creating fixed grid: %s by %s " % (mx,my)
    if abs((mx-1)*dx + xlower - xupper) > 1e-6:
        print "Warning: abs((mx-1)*dx + xlower - xupper) = ", \
              abs((mx-1)*dx + xlower - xupper)
    if abs((my-1)*dy + ylower - yupper) > 1e-6:
        print "Warning: abs((my-1)*dy + ylower - yupper) = ", \
              abs((my-1)*dy + ylower - yupper)
    tstart_max =   0.      # when to start monitoring max values
    tend_max = 1.e10       # when to stop monitoring max values
    dt_for_max = 60.       # target time increment between updating max values
    minlevel_for_max = 2   # which levels to monitor max on
    minlevel_for_arrival = 2   # which levels to monitor arrival of tsunami
    arrival_tol = 1.e-2        # tolerance for flagging arrival


    x = linspace(xlower, xupper, mx)
    y = linspace(ylower, yupper, my)
    npts = mx*my

    fname = 'fgmax_grid.txt'
    fid = open(fname,'w')
    fid.write("%g  %g             # tstart_max, tend_max\n" \
               % (tstart_max, tend_max))
    fid.write("%g                 # dt_for_max\n" % dt_for_max)
    fid.write("%i                 # minlevel_for_max\n" % minlevel_for_max)
    fid.write("%i                 # minlevel_for_arrival\n" % minlevel_for_arrival)
    fid.write("%g                 # arrival_tol\n" % arrival_tol)

    fid.write("%g  %g  %g                  # npts, mx, my \n" % (npts,mx,my))

    for j in range(my):
        for i in range(mx):
            fid.write("%20.10e %20.10e\n" % (x[i],y[j]))


    print "Created file ", fname
    fid.close()

def make_fgmax_transect():
    x1 = -120
    x2 = -70
    y1 = -17.975
    y2 = -17.975
    npts = 500
    xpts = linspace(x1,x2,npts)
    ypts = linspace(y1,y2,npts)

    print "Creating fixed grid transect with  %s points" % npts
    tstart_max =   0.      # when to start monitoring max values
    tend_max = 1.e10       # when to stop monitoring max values
    dt_for_max = 60.       # target time increment between updating max values
    minlevel_for_max = 2   # which levels to monitor max on
    minlevel_for_arrival = 2   # which levels to monitor arrival of tsunami
    arrival_tol = 1.e-2        # tolerance for flagging arrival

    fname = 'fgmax_transect.txt'
    fid = open(fname,'w')
    fid.write("%g  %g             # tstart_max, tend_max\n" \
               % (tstart_max, tend_max))
    fid.write("%g                 # dt_for_max\n" % dt_for_max)
    fid.write("%i                 # minlevel_for_max\n" % minlevel_for_max)
    fid.write("%i                 # minlevel_for_arrival\n" % minlevel_for_arrival)
    fid.write("%g                 # arrival_tol\n" % arrival_tol)

    fid.write("%g                   # npts\n" % npts)

    for j in range(npts):
        fid.write("%20.10e %20.10e\n" % (xpts[j],ypts[j]))

    print "Created file ", fname
    fid.close()

if __name__ == "__main__":
    make_fgmax_grid()
    make_fgmax_transect()


