"""
Functions to process and plot fgmax output from GeoClaw runs.

"""


from pylab import *
from clawpack.geoclaw import topotools
from numpy import ma
import os, glob, re
from matplotlib import image


#LongBeachBermGE = image.imread(WAdir + '/maps/LongBeachBermGE.png')

extent = (0,1.5,0,1.5)
x0=extent[0]
y0=extent[2]

todo_all = {
    'make h_eta.npy':  True,
    'plot h_eta':      True,
    'plot eta':        True,
    'plot speed':      True,
    'plot hmax times': True,
    'plot arrival times': True,
    'make html':       True,
    }

todo_h_eta = {
    'make h_eta.npy':  False,
    'plot h_eta':      True,
    'plot eta':        False,
    'plot speed':      False,
    'plot hmax times': False,
    'plot arrival times': False,
    'make html':       False,
    }



def make_h_eta(dirname, todo_list=todo_all, clines_t=None):
    
    outdir = dirname
    plotdir = '_plots'
    if not os.path.isdir(plotdir):
        os.system('mkdir %s' % plotdir)
        print "Created ",plotdir

    h_eta_txt = plotdir + '/h_eta.txt'
    h_eta_npy = plotdir + '/h_eta.npy'
    ert_title = plotdir

    if not os.path.isdir(dirname):
        raise Exception("Missing directory: %s" % dirname)

    if not os.path.isdir(outdir):
        raise Exception("Missing directory: %s" % outdir)

    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)


    #fname = outdir + '/fgmax_grid.txt'
    fname = 'fgmax_grid.txt'
    try:
        fid = open(fname)
    except:
        raise Exception("cannot open %s" % fname)

    # skip some lines:
    for i in range(3):
        line = fid.readline()

    line = fid.readline().split()
    fid.close()
    mx = int(line[1])
    my = int(line[2])

    fname = outdir + '/fort.FG1.valuemax' 
    print "Reading %s ..." % fname
    try:
        d = loadtxt(fname)
    except:
        raise Exception("*** Cannot read file: %s" % fname)

    h = reshape(d[:,3],(mx,my),order='F')

    fname = outdir + '/fgmax_xyB.txt'
    print "Reading %s ..." % fname
    xyB = loadtxt(fname)
    x = reshape(xyB[:,0],(mx,my),order='F')
    y = reshape(xyB[:,1],(mx,my),order='F')
    B = reshape(xyB[:,2], (mx,my),order='F')


    #import pdb; pdb.set_trace()
    eta = h + B
    h_or_eta = where(B>0, h, eta)

    thmax = reshape(d[:,7],(mx,my),order='F')  # Time maximum h recorded
    inundated = logical_and((B>0), (h>0))
    last_time_hmax = where(inundated,  thmax, 0).max()
    #last_time_hmax_hours = last_time_hmax / 3600.
    try:
        print "Latest time new maximum h recorded:  ", last_time_hmax
        #print "Latest time new maximum h recorded: %6i (%4.1f hours)" \
        #        % (last_time_hmax, last_time_hmax_hours)
    except:
        print "*** All masked? last_time_hmax = ",last_time_hmax


    # Output results as binary .npy file:

    if todo_list['make h_eta.npy']: 
        np.save(h_eta_npy, reshape(h_or_eta, (mx*my,),order='F'))
        print "Created ",h_eta_npy
    

    if todo_list['plot h_eta']:

        # Plot h or eta along with contours of topo:
        figure(1,(15,9))
        clf()
        hmax = ma.masked_where(h_or_eta==0.,h_or_eta)
        #cmax = max(hmax.max(), 20.)
        cmax = hmax.max()
        cmin = hmax.min()
        #clines = [1e-3] + list(linspace(0.5,4.5,9)) + [cmax]
        #clines = [1e-3] + list(linspace(1.0,10.0,10)) + [cmax]
        clines = linspace(cmin,cmax,10)
        colors = discrete_cmap(clines)
        contourf(x,y,hmax,clines,colors=colors)

        cbar = colorbar()
        cbar.set_ticks(clines)
        cbar.set_label('meters', fontsize=15)

        if 0:
            cbar = colorbar()
            meters = clines
            feet_labels = ['%4.1f' % (m*3.28) for m in meters]
            cbar.set_ticks(meters)
            cbar.set_ticklabels(feet_labels)
            cbar.set_label('feet', fontsize=15)


        # Contours of topo:
        clines2 = linspace(0,20,11)
        contour(x,y,B,clines2,colors='k')
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        axis(extent)
        a = gca()
        #a.set_aspect(1./cos(y0*pi/180.))
        #title(ert_title, fontsize=10)
        title("Zeta Maximum",fontsize=20)
        fname = plotdir + '/h_eta.png' 
        savefig(fname)
        print "Created ",fname


    if todo_list['plot speed']:

        # Plot max speed on top of zoom map image:
        figure(8,(13,9))
        clf()
        #imshow(LongBeachBermGE,extent=extent)
        speed = reshape(d[:,4],(mx,my),order='F')
        #speed = ma.masked_where(h_or_eta==0., speed) * 100. # convert to cm/sec
        speed = ma.masked_where(h_or_eta==0., speed) # leave as m/sec
        smax = speed.max()
        #smaxn = int(ceil(smax/100.))
        #cmax = smaxn * 100.
        #clines_s = list(linspace(0,cmax,smaxn+1))
        #cmax1 = 10.
        #if smax > 10.:
            #cmax1 = 20.
        #cmax1 = 20.
        #cmax = max(smax, cmax1)
        cmax = smax
        cmin = speed.min()
        #clines_s = list(arange(0,cmax1,cmax1/10.)) + [cmax]
        #clines_s = list(linspace(0,3,7)) + range(4,15)
        #contourf(x,y,speed,clines_s,colors=colors,alpha=0.8)
        #clines_s = [0., 1e-4, 1e-3, 1e-2, 1]
        clines_s = linspace(cmin,cmax,10)
        colors = discrete_cmap(clines_s)
        contourf(x,y,speed,clines_s)
        
        cbar = colorbar()
        cbar.set_ticks(clines_s)
        cbar.set_label('meters per second',fontsize=15)
        contour(x,y,B,[0],colors='w')
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        axis(extent)
        plot([x0],[y0],'bo')
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        axis(extent)
        a = gca()
        #a.set_aspect(1./cos(y0*pi/180.))
        title('Maximum flow speed', fontsize=20)
        fname = plotdir + '/speed_map7.png' 
        savefig(fname)
        print "Created ",fname

    if todo_list['plot hmax times']:

        # Plot time max h recorded:
        figure(9,(13,10))
        clf()
        #
        # on map:
        #imshow(LongBeachBermGE,extent=extent)
        #
        times = reshape(d[:,7],(mx,my),order='F')
        times = ma.masked_where(times < -1e50, times)      #will see the topo if masked
        if clines_t is None:
            cmax = int(ceil(times.max()))
            cmin = int(floor(times.min()))
            clines_t = list(linspace(cmin,cmax,2*(cmax-cmin)+1))
        clines_t = linspace(times.min(), times.max(), 10)
            
        times = where(h_or_eta == 0., -9999., times) # negative plots as white
        times = where(B>0, times, -9999.)                  #negative plots as white
        #times = ma.masked_where(h_or_eta == 0., times) / 3600.  # hours 
        #times = ma.masked_where(h_or_eta == 0., times) / 60.  # minutes 

        colors = discrete_cmap_times(clines_t)
        contourf(x,y,times,clines_t,colors=colors)
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('seconds',fontsize=15)

        if 0:
            # Contours of topo:
            clines2 = linspace(0,20,11)
            contour(x,y,B,clines2,colors='k')

        # on map:
        #imshow(LongBeachBermGE,extent=extent)

        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        axis(extent)
        a = gca()
        #a.set_aspect(1./cos(y0*pi/180.))
        #title(ert_title+': Time of max zeta', fontsize=20)
        plot([x0],[y0],'bo')
        title('Time of max zeta', fontsize=20)
        
        fname = plotdir + '/hmaxtimes.png' 
        savefig(fname)
        print "Created ",fname


    if todo_list['plot arrival times']:

        # Plot time max h recorded:
        figure(10,(13,9))
        clf()
        #imshow(LongBeachBermGE,extent=extent)
        atimes = reshape(d[:,11],(mx,my),order='F')
        atimes = ma.masked_where(atimes < -1e50, atimes)     #will see the topo if masked
        atimes = where(h_or_eta == 0., -9999., atimes) # negative plots as white
        atimes = where(B>0, atimes, -9999.)                  #negative plots as white
        #atimes = ma.masked_where(h_or_eta == 0., atimes) / 60.  # minutes
        #atimes = where(atimes < -1e50, -9999.,atimes)  #temporary try
        cmax = int(ceil(atimes.max()))
        cmin = int(floor(atimes.min()))
        clines_t = list(linspace(cmin,cmax,2*(cmax-cmin)+1)) 
        colors = discrete_cmap_times(clines_t)
        contourf(x,y,atimes,clines_t,colors=colors)
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('seconds',fontsize=15)

        if 0:
            # Contours of topo:
            clines2 = linspace(0,20,11)
            contour(x,y,B,clines2,colors='k')
        #
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        axis(extent)
        a = gca()
        #a.set_aspect(1./cos(y0*pi/180.))
        plot([x0],[y0],'bo')
        title('Arrival time', fontsize=20)
        fname = plotdir + '/arrival_times.png' 
        savefig(fname)
        print "Created ",fname


    return (x,y,B,hmax)

def discrete_cmap(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines. Colors go from turqouise through yellow to red.
    """
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
    Red = hstack([linspace(0,0.8,n1), ones(n2)])
    Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
    colors = zip(Red,Green,Blue)
    return colors

def discrete_cmap_times(clines):
    """
    Construct a discrete color map for the regions between the contour lines
    given in clines. For arrival times, colors go from red to turquoise.
    """
    nlines = len(clines)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = flipud(hstack([linspace(1,1,n1),linspace(1,0,n2)]))
    Red = flipud(hstack([linspace(0,0.8,n1), ones(n2)]))
    Blue = flipud(hstack([linspace(1,0.2,n1), zeros(n2)]))
    colors = zip(Red,Green,Blue)
    return colors

