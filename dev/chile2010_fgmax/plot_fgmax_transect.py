"""
Plot fgmax output from GeoClaw runs, assuming points are on transect.

"""


from pylab import *
from numpy import ma
import os

def make_plots():

    # Some things that might need to change...
    outdir = '_output'
    plotdir = '.'    # '_plots'
    fgmax_input_file = 'fgmax_transect.txt'
    #clines_zeta = None  # can set to desired contours of zeta 
    clines_zeta = [0.01] + list(linspace(0.05,0.3,6)) + [0.5,1.0,10.0]
    #clines_t = None  # can set to desired contours of arrival time or zeta time
    clines_t = linspace(0,8,17)  # hours
    clines_t_label = clines_t[::2]  # which ones to label 
    clines_t_colors = [.5,.5,.5]
    clines_topo = [0]

    plot_zeta = True
    plot_zeta_times = False
    plot_arrival_times = False



    if not os.path.isdir(outdir):
        raise Exception("Missing directory: %s" % outdir)

    if not os.path.isdir(plotdir):
        os.mkdir(plotdir)

    print outdir
    print fgmax_input_file

    # read mx and my from the input file:
    try:
        fid = open(fgmax_input_file)
    except:
        raise Exception("cannot open %s" % fgmax_input_file)

    # skip some lines:
    for i in range(5):
        line = fid.readline()

    line = fid.readline().split()
    fid.close()
    npts = int(line[0])

    fname = outdir + '/fort.FG2.valuemax' 
    print "Reading %s ..." % fname
    try:
        d = loadtxt(fname)
    except:
        raise Exception("*** Cannot read file: %s" % fname)

    x = d[:,0]
    y = d[:,1]
    y0 = 0.5*(y.min() + y.max())   # mid-latitude for scaling plots
    eta_tilde = d[:,3]

    # AMR level used for each zeta value:
    level = array(d[:,2], int)
    
    # Determine topo B at each point from the same level of AMR:
    fname = outdir + '/fort.FG2.aux1' 
    print "Reading %s ..." % fname
    daux = loadtxt(fname)
    topo = []
    for i in range(2,9):
        topoi = daux[:,i]
        topoi = ma.masked_where(topoi < -1e50, topoi)
        topo.append(topoi)

    B = ma.masked_where(level==0, topo[0])  # level==0 ==> never updated
    levelmax = level.max()
    for i in range(levelmax):
        B = where(level==i+1, topo[i], B)

    h = where(eta_tilde > B, eta_tilde - B, 0.)


    # zeta = max h on land or max eta offshore:
    zeta = where(B>0, h, eta_tilde)

    tzeta = d[:,7]  # Time maximum h recorded
    tzeta = ma.masked_where(tzeta < -1e50, tzeta)      
    tzeta = ma.masked_where(zeta == 0., tzeta) / 3600.  # hours 

    inundated = logical_and((B>0), (h>0))

    atimes = d[:,11]
    atimes = ma.masked_where(atimes < -1e50, atimes)  
    atimes = ma.masked_where(zeta == 0., atimes) / 3600.  # hours 

    if plot_zeta:

        # Plot h or eta along with contours of topo:
        figure(101)
        clf()
        zeta = ma.masked_where(zeta==0.,zeta)
        subplot(211)
        plot(x,zeta)
        title("Zeta Maximum at y = %g" % y[0])

        subplot(212)
        plot(x,zeta)
        plot(x,B,'g')  # topo
        title("Topography")
        
        fname = plotdir + '/zeta_transect.png' 
        savefig(fname)
        print "Created ",fname


    if plot_zeta_times:

        # Plot time max h recorded:
        figure(102)
        clf()

        if clines_t is None:
            clines_t = linspace(tzeta.min(), tzeta.max(), 10)
            

        colors = discrete_cmap_times(clines_t)
        contourf(x,y,tzeta,clines_t,colors=colors)
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('hours',fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='k',linestyles='-')

        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))
        title('Time of max zeta', fontsize=20)
        
        fname = plotdir + '/zetatimes.png' 
        savefig(fname)
        print "Created ",fname


    if plot_arrival_times:

        # Plot time max h recorded:
        figure(103)
        clf()
        if clines_t is None:
            clines_t = linspace(atimes.min(), atimes.max(), 10)
        colors = discrete_cmap_times(clines_t)
        contourf(x,y,atimes,clines_t,colors=colors)
        cbar = colorbar()
        cbar.set_ticks(clines_t)
        cbar.set_label('hours',fontsize=15)

        # Contours of topo:
        contour(x,y,B,clines_topo,colors='k',linestyles='-')
        #
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        gca().set_aspect(1./cos(y0*pi/180.))
        title('Arrival time', fontsize=20)
        fname = plotdir + '/arrival_times.png' 
        savefig(fname)
        print "Created ",fname



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

if __name__ == "__main__":
    make_plots()
