from clawpack.geoclaw import kmltools
import os
from numpy import sqrt

class fgmax_grid_parameters(object):

    def __init__(self):
        self.point_style = 2      # expect 1 or 2
        self.x1 = None
        self.x2 = None
        self.y1 = None
        self.y2 = None
        self.dx = None
        self.dy = None
        self.nx = None
        self.ny = None
        self.tstart_max =  0.
        self.tend_max = 1.e10     # when to stop monitoring max values
        self.dt_check = 10.       # target time (sec) increment between updating 
                                  # max values
        self.min_level_check = None    # which levels to monitor max on
        self.arrival_tol = 1.e-2    # tolerance for flagging arrival
        self.fname = 'fgmax_grid.txt'

def make_fgmax(FG):
    print "---------------------------------------------- "
    x1,x2 = FG.x1, FG.x2
    y1,y2 = FG.y1, FG.y2
    point_style = FG.point_style
    if point_style not in [1,2]:
        raise NotImplementedError("make_fgmax not implemented for point_style %i" \
            % point_style)

    if point_style == 2:
        # 2d grid of points
        if FG.nx is None:
            dx = FG.dx
            nx = int(round((x2-x1)/dx)) + 1  
            if abs((nx-1)*dx + x1 - x2) > 1e-6:
                print "Warning: abs((nx-1)*dx + x1 - x2) = ", \
                      abs((nx-1)*dx + x1 - x2)
                print "         old x2: %22.16e" % x2
                x2 = x1 + dx*(nx-1)
                print "         resetting x2 to %22.16e" % x2
        else:
            nx = FG.nx
            dx = (x2-x1)/(nx+1.)
            if FG.dx is not None:
                print "*** Warning: dx specified over-ridden by: ",dx
    
        if FG.ny is None:
            dy = FG.dy
            if dy is None:
                dy = dx
            ny = int(round((y2-y1)/dy)) + 1  
            if abs((ny-1)*dy + y1 - y2) > 1e-6:
                print "Warning: abs((ny-1)*dy + y1 - y2) = ", \
                      abs((ny-1)*dy + y1 - y2)
                print "         old y2: %22.16e" % y2
                y2 = y1 + dy*(ny-1)
                print "         resetting y2 to %22.16e" % y2
        else:
            ny = FG.ny
            dy = (y2-y1)/(ny+1.)
            if FG.dy is not None:
                print "*** Warning: dy specified over-ridden by: ",dy
    
    
        npts = nx*ny
    
        fid = open(FG.fname,'w')
        fid.write("%16.10e            # tstart_max\n"  % FG.tstart_max)
        fid.write("%16.10e            # tend_max\n"  % FG.tend_max)
        fid.write("%16.10e            # dt_check\n" % FG.dt_check)
        fid.write("%i %s              # min_level_check\n" \
                            % (FG.min_level_check,16*" "))

        fid.write("%16.10e            # arrival_tol\n" % FG.arrival_tol)
        fid.write("%i %s              # point_style\n" \
                            % (FG.point_style,16*" "))
    
        fid.write("%i  %i %s          # nx,ny\n" \
                            % (nx,ny,10*" "))
        fid.write("%16.10e   %20.10e            # x1, y1\n" % (x1,y1))
        fid.write("%16.10e   %20.10e            # x2, y2\n" % (x2,y2))
        fid.close()
        
    
        print "Created file ", FG.fname
        print "   specifying fixed grid with shape %i by %i, with  %i points" \
                % (nx,ny,npts)
        print "   lower left  = (%15.10f,%15.10f)" % (x1,y1)
        print "   upper right = (%15.10f,%15.10f)" % (x2,y2)
        print "   dx = %15.10e,  dy = %15.10e" % (dx,dy)
    
        xy = [x1,x2,y1,y2]
        fname_root = os.path.splitext(FG.fname)[0]
        kml_file = fname_root + '.kml'
        kmltools.box2kml(xy, kml_file, fname_root, color='8888FF')

    elif point_style==1:
        # 1d transect of points
        if FG.npts is None:
            dx = FG.dx
            npts = int(round(sqrt((x2-x1)**2 + (y2-y1)**2)/dx)) + 1
            if abs((npts-1)*dx + x1 - x2) > 1e-6:
                print "Warning: abs((npts-1)*dx + x1 - x2) = ", \
                      abs((npts-1)*dx + x1 - x2)
                x2 = x1 + dx*(npts-1)
                y2 = y1 + dx*(npts-1)
                print "         resetting x2 to %g" % x2
                print "         resetting y2 to %g" % y2
        else:
            npts = FG.npts
            dx = sqrt((x2-x1)**2 + (y2-y1)**2)/(npts+1.)
            if FG.dx is not None:
                print "*** Warning: dx specified over-ridden by: ",dx
    
    
        print "Creating 1d fixed grid with %s points" % npts
        print "   dx = %g" % dx
    
        fid = open(FG.fname,'w')
        fid.write("%g                 # tstart_max\n"  % FG.tstart_max)
        fid.write("%g                 # tend_max\n"  % FG.tend_max)
        fid.write("%g                 # dt_check\n" % FG.dt_check)
        fid.write("%i                 # min_level_check\n" % FG.min_level_check)
        fid.write("%g                 # arrival_tol\n" % FG.arrival_tol)
        fid.write("%g                 # point_style\n" % FG.point_style)
    
        fid.write("%i                 # npts\n" % (npts))
        fid.write("%g   %g            # x1, y1\n" % (x1,y1))
        fid.write("%g   %g            # x2, y2\n" % (x2,y2))
        fid.close()
        
    
        print "Created file ", FG.fname
        print "   specifying fixed grid with %i points equally spaced from " \
                % npts
        print "   (%g,%g)  to  (%g,%g)" % (x1,y1,x2,y2)
        
        # not yet implemented:
        #fname_root = os.path.splitext(FG.fname)[0]
        #kml_file = fname_root + '.kml'
        #kmltools.line2kml(xy, kml_file, fname_root, color='8888FF')

