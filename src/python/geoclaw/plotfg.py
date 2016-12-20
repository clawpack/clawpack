"""
Plot fixed grid output.
fgno refers to the number of the fixed grid

Now creates html pages of plots by default if
executed from Unix command line via:
    $ python plotfg.py 1
to plot fixed grid number 1, for example.

Specify output directory other than _output by adding the outdir to this line:
    $ python plotfg.py 1 _output-05

"""

from __future__ import absolute_import
from __future__ import print_function
from pylab import *
from clawpack.visclaw import geoplot, colormaps, plotpages
# Cannot find this!
# from plottools import fix_long_tick_labels

import os
from numpy import ma
from clawpack.clawutil.data import ClawData
from six.moves import range
from six.moves import input


class ClawPlotFGData(ClawData):

    def __init__(self,fgno=1):
        super(ClawPlotFGData,self).__init__()
        self.add_attribute('fgno',fgno)
        self.add_attribute('drytol',1.e-2)
        self.add_attribute('exposed_tol',1.e-2)
        self.add_attribute('clines',linspace(-20,20,21))

        self.add_attribute('water_cmap', geoplot.tsunami_colormap)
        self.add_attribute('land_cmap', geoplot.land_colors)
        self.add_attribute('seafloor_cmap', geoplot.seafloor_colormap)
        self.add_attribute('inundated_cmap', \
             colormaps.make_colormap({0.:[1,.9,.9], 1.:[1,0,0]}))

        self.add_attribute('water_clim',(-1,1))
        self.add_attribute('land_clim',(0,10))
        self.add_attribute('seafloor_clim',(-1,0))
        self.add_attribute('inundated_clim',(0,1))

        self.add_attribute('water_add_colorbar',True)
        self.add_attribute('land_add_colorbar',False)
        self.add_attribute('seafloor_add_colorbar',True)
        self.add_attribute('inundated_add_colorbar',True)

        self.add_attribute('eta_show',True)
        self.add_attribute('seafloor_show',True)
        self.add_attribute('inundated_show',True)

        self.add_attribute('outdir','_output')
        self.add_attribute('plotdir',None)
        self.add_attribute('save_png',False)
        self.add_attribute('solutions',{})
        self.add_attribute('grids',{})
        self.add_attribute('combined_figure',True)

    def list_frames(self):

        import glob
        pattern = "%s/fort.fg%s_*" % (self.outdir,str(self.fgno).zfill(2))
        files = glob.glob(pattern)
        if len(files) == 0:
            print('*** No files found of form ', pattern)
        framenos = []
        for file in files:
            line = open(file,'r').readline()
            t = float(line.split()[0])
            print("%s: t = %s" % (file,t))
            frameno = file[-2:]
            framenos.append(int(frameno))
        return framenos


    def get_frame(self, frameno):

        if frameno in self.solutions:
            # don't read if already in dictionary:
            return self.grid, self.solutions[frameno]

        fname = "fort.fg%s_%s" % (str(self.fgno).zfill(2), str(frameno).zfill(4))
        fname = os.path.join(self.outdir,fname)
        if not os.path.exists(fname):
            print("*** Did not find file ",fname," in directory ",self.outdir)
            raise IOError("Missing fixed grid output file")
        

        self.grid = ClawData()
        grid = self.grid

        # Read parameters from header:

        file = open(fname,'r')

        line = file.readline()
        t = float(line.split()[0])
        print("Reading fixed grid output from ",fname)
        print('   Frame %s at t = %s' % (frameno,t))

        line = file.readline()
        grid.mx = int(line.split()[0])

        line = file.readline()
        grid.my = int(line.split()[0])

        line = file.readline()
        grid.xlow = float(line.split()[0])

        line = file.readline()
        grid.ylow = float(line.split()[0])

        line = file.readline()
        grid.xhi = float(line.split()[0])

        line = file.readline()
        grid.yhi = float(line.split()[0])

        file.close()

        grid.x = linspace(grid.xlow, grid.xhi, grid.mx+1)
        grid.y = linspace(grid.ylow, grid.yhi, grid.my+1)
        grid.dx = grid.x[1]-grid.x[0]
        grid.dy = grid.y[1]-grid.y[0]
        grid.xcenter = grid.x[:-1] + grid.dx/2.
        grid.ycenter = grid.y[:-1] + grid.dy/2.

        d = loadtxt(fname, skiprows=8)

        solution = ClawData()
        solution.t = t
        solution.ncols = d.shape[1]
        solution.h = reshape(d[:,0], (grid.my,grid.mx))
        solution.B = reshape(d[:,3], (grid.my,grid.mx))
        solution.eta = reshape(d[:,4], (grid.my,grid.mx))
        solution.surface = ma.masked_where(isnan(solution.eta),solution.eta)
        solution.land = ma.masked_where(solution.h>self.drytol,solution.B)
        solution.fg = empty((solution.ncols,grid.my,grid.mx), dtype=float)
        for col in range(solution.ncols):
            solution.fg[col,:,:] = reshape(d[:,col],(grid.my,grid.mx))
        
        self.solutions[frameno] = solution
        return grid, solution


    def plotfg(self, frameno):

        grid, solution = self.get_frame(frameno)
        print("Plotting frame %s at time t = %s"  % (frameno,solution.t))

    
        # Define function to plot topo contours for use in multiple places:
        def add_contours():
            contour(grid.xcenter,grid.ycenter,solution.B,self.clines,colors='k')
            contour(grid.xcenter,grid.ycenter,solution.B,[0.],colors='k',linewidths=2)
    
        #t_str = "%10.3f" % solution.t
        tmin = solution.t / 60.
        t_str = "%8.2f minutes" % tmin
    
        if self.combined_figure:
            figno = 153
            figure(figno, (14,8))
            clf()
            subplot(121)
        else:
            figno = 150
            figure(figno)
            clf()
        
        if ma.count(solution.surface) != 0:
            pcolormesh(grid.x,grid.y,solution.surface,cmap=self.water_cmap)
            clim(self.water_clim)
            if self.water_add_colorbar: colorbar(shrink=0.5)
        
        if ma.count(solution.land) != 0:
            pcolormesh(grid.x,grid.y,solution.land,cmap=self.land_cmap)
            clim(self.land_clim)
            if self.land_add_colorbar: colorbar(shrink=0.5)
        
        add_contours()
        
        
        title('Eta on FG %s at time = %s' % (self.fgno,t_str))
        xlim((grid.xlow, grid.xhi))
        ylim((grid.ylow, grid.yhi,))
        # fix_long_tick_labels()        
        axis('tight')
        axis('scaled')
    
        if self.save_png:
            fname = 'FixedGrid%sFrame%sfig%s.png' \
                %  (str(self.fgno).zfill(2), str(frameno).zfill(4), figno)
            savefig(fname)
            print("Saved figure as ",fname)
        
        
        if solution.ncols > 5:

            etamin = solution.fg[5,:,:]
            etamax = solution.fg[6,:,:]
    
            etamax2 = where(solution.B<0, 1., etamax)
            
            # Add red contour of maximum eta
            #contour(xcenter,ycenter,etamax2,[drytol],colors='r',linewidths=2)
            
            # Add brown contour of minimum eta
            #contour(xcenter,ycenter,etamin-B,[exposed_tol],colors=([.9,.8,.2],),linewidths=2)
    
            # Determine exposed and inundatated regions:
    
            B = solution.B
            exposed_tol = self.exposed_tol
            drytol = self.drytol
            exposed = ma.masked_where(((B>0) | (etamin > B+exposed_tol)), etamin)
            inundated = ma.masked_where(((B<0) | (etamax < B+drytol)), etamax-B)

        #---------------------------------------------------------------

        if self.combined_figure:
            subplot(122)

        if self.inundated_show:
            if solution.ncols < 7:
                raise ValueError("*** Data does not include etamin/etamax")
    
            if not self.combined_figure:
                figno = 151
                figure(figno)
                clf()
            
            if ma.count(inundated) != 0:
                pcolormesh(grid.x,grid.y,inundated,cmap=self.inundated_cmap)
                clim(self.inundated_clim)
                if self.inundated_add_colorbar: colorbar(shrink=0.5)
            add_contours()
            # Add red contour of maximum eta
            #contour(xcenter,ycenter,etamax2,[drytol],colors='r',linewidths=2)
            title("Inundated region for t <= %s" % t_str)
            xlim((grid.xlow, grid.xhi))
            ylim((grid.ylow, grid.yhi,))
            # fix_long_tick_labels()        
            axis('tight')
            axis('scaled')
            if self.save_png:
                fname = 'FixedGrid%sFrame%sfig%s.png' \
                    %  (str(self.fgno).zfill(2), str(frameno).zfill(4), figno)
                savefig(fname)
                print("Saved figure as ",fname)
        
        #---------------------------------------------------------------

        if self.seafloor_show:
            if solution.ncols < 7:
                raise ValueError("*** Data does not include etamin/etamax")
    
            if not self.combined_figure:
                figno = 152
                figure(figno)
                clf()
            
            if ma.count(exposed) != 0:
                pcolormesh(grid.x,grid.y,exposed,cmap=self.seafloor_cmap)
                clim(self.seafloor_clim)
                if self.seafloor_add_colorbar: colorbar(shrink=0.5)
            add_contours()
            # Add brown contour of minimum eta
            #contour(xcenter,ycenter,etamin-B,[exposed_tol],colors=([.9,.8,.2],),linewidths=2)
            if self.combined_figure:
                title("Eta min/max t <= %s" % t_str)
            else:
                title("Exposed seabed for t <= %s" % t_str)
            xlim((grid.xlow, grid.xhi))
            ylim((grid.ylow, grid.yhi,))
            # fix_long_tick_labels()        
            axis('tight')
            axis('scaled')
            if self.save_png:
                fname = 'FixedGrid%sFrame%sfig%s.png' \
                    %  (str(self.fgno).zfill(2), str(frameno).zfill(4), figno)
                savefig(fname)
                print("Saved figure as ",fname)
                
            
        
    def fgloop(self):
        for frameno in range(1,100):
            try:
                self.plotfg(frameno)
            except IOError:
                break
            ans = input("Hit return for next time, q to quit, s to savefig... ")
            if ans=='s':
                fname = 'FixedGrid%sFrame%s.png' %  (str(self.fgno).zfill(2), str(frameno).zfill(4))
                savefig(fname)
                print("Saved figure as ",fname)
                ans = input("Hit return for next time, q to quit, s to savefig... ")
                
            if ans=='q':
                break
    
                
    def fg2html(self,framenos='all'):
        if self.plotdir is None:
            self.plotdir='_fgplots_fg%s' % str(self.fgno).zfill(2)
        startdir = os.getcwd()
        self.outdir = os.path.abspath(self.outdir)
        plotpages.cd_with_mkdir(self.plotdir, overwrite=True)
        self.save_png = True

        ppd = plotpages.PlotPagesData()
        ppd.timeframes_frametimes = {}

        if framenos=='all':
            framenos = list(range(1,200))
    
        for frameno in framenos:
            try:
                grid, solution = self.get_frame(frameno)
                t = solution.t
                ppd.timeframes_frametimes[frameno] = t
                self.plotfg(frameno)
            except IOError:
                break
            draw()
    
        os.chdir(startdir)
        ppd.plotdir = self.plotdir
        ppd.html_index_fname = "_PlotIndex_FixedGrid%s.html" \
            % str(self.fgno).zfill(2)
        ppd.html_index_title = "Fixed Grids Plot Index"
        ppd.timeframes_prefix='FixedGrid%sFrame' % str(self.fgno).zfill(2)
        ppd.html_homelink = "_PlotIndex.html"
        ppd.timeframes_fignames[150] = 'Surface'
        ppd.timeframes_fignames[151] = 'Inundation'
        ppd.timeframes_fignames[152] = 'Exposed seafloor'
        plotpages.timeframes2html(ppd)
            


if __name__ == "__main__":
    # if executed at the Unix command line....
    import sys
    args = sys.argv[1:]   # any command line arguments
    fg = ClawFGData()
    if len(args) > 0:
        fg.fgno = args[0]
    if len(args) > 1:
        fg.outdir = args[1]
    fg.fg2html(fg)
