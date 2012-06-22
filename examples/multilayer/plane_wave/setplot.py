
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

def meters_to_km(x,y):
    return x/1e3,y/1e3
    
#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    import os

    import numpy as np
    import matplotlib.pyplot as plt

    from visclaw.plotters import colormaps, geoplot
    from oldclawdata import Data

    amrdata = Data(os.path.join(plotdata.outdir,'amr2ez.data'))
    hurricane_data = Data(os.path.join(plotdata.outdir,'hurricane.data'))
    multilayer_data = Data(os.path.join(plotdata.outdir,'multilayer.data'))
    
    def transform_c2p(x,y,x0,y0,theta):
        return ((x+x0)*np.cos(theta) - (y+y0)*np.sin(theta),
                (x+x0)*np.sin(theta) + (y+y0)*np.cos(theta))

    def transform_p2c(x,y,x0,y0,theta):
        return ( x*np.cos(theta) + y*np.sin(theta) - x0,
                -x*np.sin(theta) + y*np.cos(theta) - y0)
        
    
    if multilayer_data.bathy_type == 1:
        x = [0.0,0.0]
        y = [0.0,1.0]
        x1,y1 = transform_c2p(x[0],y[0],multilayer_data.bathy_location,
                                    0.0,multilayer_data.bathy_angle)
        x2,y2 = transform_c2p(x[1],y[1],multilayer_data.bathy_location,
                                    0.0,multilayer_data.bathy_angle)
        
        m = (y1 - y2) / (x1 - x2)
        x[0] = (amrdata.ylower - y1) / m + x1
        y[0] = amrdata.ylower
        x[1] = (amrdata.yupper - y1) / m + x1
        y[1] = amrdata.yupper
        ref_lines = [((x[0],y[0]),(x[1],y[1]))]
    else:
        ref_lines = []
    print ref_lines
    
    plotdata.clearfigures()
    plotdata.clear_frames = False
    plotdata.clear_figs = True
    
    plotdata.save_frames = False
    
    # ========================================================================
    #  Generic helper functions
    # ========================================================================
    def pcolor_afteraxes(current_data):
        bathy_ref_lines(current_data)
        gauge_locations(current_data)
        
    def contour_afteraxes(current_data):
        # gauge_locations(current_data)
        m_to_km_labels()
        # plt.hold(True)
        # pos = -80.0 * (23e3 / 180) + 500e3 - 5e3
        # plt.plot([pos,pos],[-300e3,300e3],'b',[pos-5e3,pos-5e3],[-300e3,300e3],'y')
        # plt.hold(False)
        # wind_contours(current_data)
        # bathy_ref_lines(current_data)
        
    def profile_afteraxes(current_data):
        pass
        
    def bathy_ref_lines(current_data):
        plt.hold(True)
        for ref_line in ref_lines:
            x1 = ref_line[0][0]
            y1 = ref_line[0][1]
            x2 = ref_line[1][0]
            y2 =  ref_line[1][1]
            plt.plot([x1,x2],[y1,y2],'y--',linewidth=1)
        plt.hold(False)
        
    def gauge_locations(current_data,gaugenos='all'):
        from pyclaw.plotters import gaugetools
        plt.hold(True)
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=gaugenos, format_string='kx', add_labels=True)
        plt.hold(False)
    # ========================================================================
    #  Data extraction routines
    #     0    1     2     3    4     5      6     7      8      9
    #   h(1),hu(1),hv(1),h(2),hu(2),hv(2),eta(1),eta(2),wind_x,wind_y    
    # ========================================================================    

    def extract_eta(h,eta,DRY_TOL=10**-3):
        index = np.nonzero((np.abs(h) < DRY_TOL) + (h == np.nan))
        eta[index[0],index[1]] = np.nan
        return eta
        
    def eta1(cd):
        # return cd.q[:,:,6]
        return extract_eta(cd.q[:,:,0],cd.q[:,:,6])
        
    def eta2(cd):
        # return cd.q[:,:,7]
        return extract_eta(cd.q[:,:,3],cd.q[:,:,7])
    
    def b(current_data):
        h1 = current_data.q[:,:,0]
        h2 = current_data.q[:,:,3]
        
        return current_data.q[:,:,6] - h1 - h2
    
    def extract_velocity(h,hu,DRY_TOL=10**-8):
        u = np.ones(hu.shape) * np.nan
        # u = np.zeros(hu.shape)
        index = np.nonzero((np.abs(h) > DRY_TOL) * (h != np.nan))
        u[index[0],index[1]] = hu[index[0],index[1]] / h[index[0],index[1]]
        return u
    
    def water_u1(cd):
        return extract_velocity(cd.q[:,:,0],cd.q[:,:,1])
        
        # return np.where(abs(current_data.q[:,:,0]) > 10**-16,
        #     current_data.q[:,:,1] / current_data.q[:,:,0],0.0)
            
    def water_u2(cd):
        return extract_velocity(cd.q[:,:,3],cd.q[:,:,4])
        # return np.where(abs(current_data.q[:,:,3]) > 10**-16,
        #     current_data.q[:,:,4] / current_data.q[:,:,3],0.0)
        
    def water_v1(cd):
        return extract_velocity(cd.q[:,:,0],cd.q[:,:,2])
        # return np.where(abs(current_data.q[:,:,0]) > 10**-16,
        #     current_data.q[:,:,2] / current_data.q[:,:,0],0.0)
            
    def water_v2(cd):
        return extract_velocity(cd.q[:,:,3],cd.q[:,:,5])
        # return np.where(abs(current_data.q[:,:,3]) > 10**-16,
        #     current_data.q[:,:,5] / current_data.q[:,:,3],0.0)
        
    def water_speed1(current_data):
        u = water_u1(current_data)
        v = water_v1(current_data)
            
        return np.sqrt(u**2+v**2)
        
    def water_speed2(current_data):
        u = water_u2(current_data)
        v = water_v2(current_data)
            
        return np.sqrt(u**2+v**2)
        
    def water_speed_depth_ave(current_data):
        h1 = current_data.q[:,:,0]
        h2 = current_data.q[:,:,3]
        u1 = water_speed1(current_data)
        u2 = water_speed1(current_data)
        
        return (h1*u1 + h2*u2) / (h1+h2)
        
    def water_quiver1(current_data):
        u = water_u1(current_data)
        v = water_v1(current_data)
            
        plt.hold(True)
        Q = plt.quiver(current_data.x[::2,::2],current_data.y[::2,::2],
                        u[::2,::2],v[::2,::2])
        max_speed = np.max(np.sqrt(u**2+v**2))
        label = r"%s m/s" % str(np.ceil(0.5*max_speed))
        plt.quiverkey(Q,0.15,0.95,0.5*max_speed,label,labelpos='W')
        plt.hold(False)
        
    def water_quiver2(current_data):
        u = water_u1(current_data)
        v = water_v1(current_data)
            
        plt.hold(True)
        Q = plt.quiver(current_data.x[::2,::2],current_data.y[::2,::2],
                        u[::2,::2],v[::2,::2])
        max_speed = np.max(np.sqrt(u**2+v**2))
        label = r"%s m/s" % str(np.ceil(0.5*max_speed))
        plt.quiverkey(Q,0.15,0.95,0.5*max_speed,label,labelpos='W')
        plt.hold(False)
    
    # ========================================================================
    #  Plot items
    # ========================================================================
    def add_surface_elevation(plotaxes,surface,bounds=None,plot_type='pcolor'):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            # plotitem.plot_var = geoplot.surface
            if surface == 1:
                plotitem.plot_var = eta1
            elif surface == 2:
                plotitem.plot_var = eta2
            if bounds is not None:                
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            # plotitem.pcolor_cmap = geoplot.tsunami_colormap
            plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]

        elif plot_type == 'contour':
            plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
            plotitem.plot_var = surface + 5
            if bounds is not None:
                plotitem.contour_levels = bounds
            plotitem.amr_contour_show = [1,1,1]
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
            plotitem.amr_contour_colors = 'k'
        else:
            raise NotImplementedError("Plot type %s not implemented" % plot_type)
        
    def add_layer_depth(plotaxes,layer,bounds=None,plot_type='pcolor'):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            if layer == 1:
                plotitem.plot_var = 0
            elif layer == 2:
                plotitem.plot_var = 3
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.imshow_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
    
    def add_speed(plotaxes,layer,bounds=None,plot_type='pcolor'):        
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            if layer == 1:
                plotitem.plot_var = water_speed1
            elif layer == 2:
                plotitem.plot_var = water_speed2
            plotitem.imshow_cmap = plt.get_cmap('PuBu')
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1]
        elif plot_type == 'contour':
            pass

    def add_x_velocity(plotaxes,layer,plot_type='pcolor',bounds=None):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            if layer == 1:
                plotitem.plot_var = water_u1
            if layer == 2:
                plotitem.plot_var = water_u2
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.add_colorbar = True
            plotitem.imshow_cmap = plt.get_cmap('PiYG')
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1]
        elif plot_type == 'contour':
            pass
    
    def add_y_velocity(plotaxes,layer,plot_type='pcolor',bounds=None):
        if plot_type == 'pcolor' or plot_type == 'imshow':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            if layer == 1:
                plotitem.plot_var = water_v1
            if layer == 2:
                plotitem.plot_var = water_v2
            if bounds is not None:
                plotitem.imshow_cmin = bounds[0]
                plotitem.imshow_cmax = bounds[1]
            plotitem.imshow_cmap = plt.get_cmap('PiYG')
            plotitem.add_colorbar = True
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1]
        elif plot_type == 'contour':
            pass
    
    
    # Land
    def add_land(plotaxes,plot_type='pcolor'):
        r"""Add plot item for land"""
        
        if plot_type == 'pcolor':
            plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
            plotitem.show = True
            plotitem.plot_var = geoplot.land
            plotitem.imshow_cmap = geoplot.land_colors
            plotitem.imshow_cmin = 0.0
            plotitem.imshow_cmax = 80.0
            plotitem.add_colorbar = False
            plotitem.amr_gridlines_show = [0,0,0]
            plotitem.amr_gridedges_show = [1,1,1]
        elif plot_type == 'contour':            
            plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
            plotitem.plot_var = geoplot.land
            plotitem.contour_nlevels = 40
            plotitem.contour_min = 0.0
            plotitem.contour_max = 100.0
            plotitem.amr_contour_colors = ['g']  # color on each level
            plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
            plotitem.gridlines_show = 0
            plotitem.gridedges_show = 0
            plotitem.show = True
        else:
            raise NotImplementedError("Plot type %s not implemented" % plot_type)
            
    def add_combined_profile_plot(plot_data,slice_value,direction='x',figno=120):
        def slice_index(cd):
            if direction == 'x':
                if cd.grid.y.lower < slice_value < cd.grid.y.upper:
                    return int((slice_value - cd.grid.y.lower) / cd.dy - 0.5)
                else:
                    return None
            elif direction == 'y':
                if cd.grid.x.lower < slice_value < cd.grid.x.upper:
                    return int((slice_value - cd.grid.x.lower) / cd.dx - 0.5)
                else:
                    return None
    
        def bathy_profile(current_data):
            index = slice_index(current_data)
            if direction == 'x':
                if index:
                    return current_data.x[:,index], b(current_data)[:,index]
                else:
                    return None
            elif direction == 'y':
                if index:
                    return current_data.y[index,:], b(current_data)[index,:]
                else:
                    return None
    
        def lower_surface(current_data):
            index = slice_index(current_data)
            if direction == 'x':
                if index:
                    return current_data.x[:,index], eta2(current_data)[:,index]
                else:
                    return None
            elif direction == 'y':
                if index:
                    return current_data.y[index,:], eta2(current_data)[index,:]
                else:
                    return None
        
        def upper_surface(current_data):
            index = slice_index(current_data)
            if direction == 'x':
                if index:
                    return current_data.x[:,index], eta1(current_data)[:,index]
                else:
                    return None
            elif direction == 'y':
                if index:
                    return current_data.y[index,:], eta1(current_data)[index,:]
                else:
                    return None
        
        # Surfaces
        plotfigure = plotdata.new_plotfigure(name='combined_surface_%s' % figno,figno=figno)
        plotfigure.show = True
        plotfigure.kwargs = {'figsize':(6,6)}
    
        # Top surface
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.axescmd = 'subplot(2,1,1)'
        plotaxes.title = 'Surfaces Profile %s at %s' % (direction,slice_value)
        if multilayer_data.init_type == 2:
            plotaxes.xlimits = xlimits
        elif multilayer_data.init_type == 6:
            plotaxes.xlimits = ylimits
    
        plotaxes.ylimits = top_surf_zoomed
        def top_surf_afteraxes(cd):
            plt.hold(True)
            plt.xlabel('')
            locs,labels = plt.xticks()
            labels = ['' for i in xrange(len(locs))]
            plt.xticks(locs,labels)
            plt.plot([multilayer_data.bathy_location,multilayer_data.bathy_location],top_surf_zoomed,'--k')
            plt.ylabel('m')
            plt.hold(False)
        plotaxes.afteraxes = top_surf_afteraxes
        plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
        plotitem.map_2d_to_1d = upper_surface
        plotitem.amr_plotstyle = ['-','+','x']
        # plotitem.color = (0.2,0.8,1.0)
        plotitem.show = True
    
        # Internal surface
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.axescmd = 'subplot(2,1,2)'
        plotaxes.title = ''
        if multilayer_data.init_type == 2:
            plotaxes.xlimits = xlimits
        elif multilayer_data.init_type == 6:
            plotaxes.xlimits = ylimits
        plotaxes.ylimits = bottom_surf_zoomed
        def internal_surf_afteraxes(cd):
            plt.hold(True)
            plt.title('')
            plt.ylabel('m')
            plt.subplots_adjust(hspace=0.05)        
            plt.plot([multilayer_data.bathy_location,multilayer_data.bathy_location],bottom_surf_zoomed,'--k')
            plt.hold(False)
        plotaxes.afteraxes = internal_surf_afteraxes
        plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
        plotitem.map_2d_to_1d = lower_surface
        plotitem.amr_plotstyle = ['-','+','x']
        plotitem.color = 'k'
        plotitem.show = True  
    
    def add_velocities_profile_plot(plot_data,slice_value,direction='x',figno=130):
        
        def slice_index(cd):
            if cd.grid.y.lower < slice_value < cd.grid.y.upper:
                return int((slice_value - cd.grid.y.lower) / cd.dy - 0.5)
            else:
                return None
                
        
        def upper_surface(current_data):
            index = slice_index(current_data)
            if index:
                return current_data.x[:,index], eta1(current_data)[:,index]
            else:
                return None
        
        def top_speed(current_data):
            index = slice_index(current_data)
            if index:
                return current_data.x[:,index],water_u1(current_data)[:,index]
            else:
                return None, None
                        
        def bottom_speed(current_data):
            index = slice_index(current_data)
            if index:
                return current_data.x[:,index],water_u2(current_data)[:,index]
            else:
                return None, None
        
        # Velocities
        plotfigure = plotdata.new_plotfigure(name='combined_velocities_%s' % figno,figno=figno)
        plotfigure.show = True
        # plotfigure.kwargs = {'figsize':(6,6)}
    
        # Top surface
        plotaxes = plotfigure.new_plotaxes()
        # plotaxes.axescmd = 'subplot(2,1,1)'
        plotaxes.title = 'Velocities Profile %s at %s' % (direction,slice_value)
        if multilayer_data.init_type == 2:
            plotaxes.xlimits = xlimits
        elif multilayer_data.init_type == 6:
            plotaxes.xlimits = ylimits
        plotaxes.ylimits = velocities_zoomed
        def velocity_afteraxes(cd):
            plt.hold(True)
            plt.xlabel('')
            locs,labels = plt.xticks()
            labels = ['' for i in xrange(len(locs))]
            plt.xticks(locs,labels)
            plt.plot([multilayer_data.bathy_location,multilayer_data.bathy_location],velocities_zoomed,'--k')
            plt.ylabel('m/s')
            plt.hold(False)
        plotaxes.afteraxes = velocity_afteraxes
    
        plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
        plotitem.map_2d_to_1d = top_speed
        plotitem.amr_plotstyle = ['-','+','x']
        # plotitem.color = (0.2,0.8,1.0)
        plotitem.show = True
    
        plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
        plotitem.map_2d_to_1d = bottom_speed
        plotitem.amr_plotstyle = ['-','+','x']
        plotitem.color = 'k'
        plotitem.show = True  

    # ========================================================================
    # Axis limits
    #xlimits = [amrdata.xlower,amrdata.xupper]
    xlimits = [-0.5,0.5]
    xlimits_zoomed = xlimits
    #ylimits = [amrdata.ylower,amrdata.yupper]
    ylimits = [-0.5,0.5]
    eta = [multilayer_data.eta[0],multilayer_data.eta[1]]
    top_surface_limits = [eta[0]-0.03,eta[0]+0.03]
    internal_surface_limits = [eta[1]-0.015,eta[1]+0.015]
    top_speed_limits = [0.0,0.1]
    internal_speed_limits = [0.0,0.03]
    
    top_surf_zoomed = [eta[0] - 0.03,eta[0]+0.03]
    bottom_surf_zoomed = [-0.7,-0.5]
    # bottom_surf_zoomed = [eta[1] - 0.5,eta[1] + 0.5]
    velocities_zoomed = [-0.005,0.005]
    
    # Single layer test limits
    # top_surface_limits = [eta[0]-2.5,eta[0]+2.5]
    # top_speed_limits = [0.0,6.0]
    
    # ========================================================================
    #  Surface Elevations
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface', figno=0)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize':(14,4)}
    
    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Top Surface'
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    add_surface_elevation(plotaxes,1,bounds=top_surface_limits)
    # add_surface_elevation(plotaxes,1,bounds=[-0.06,0.06])
    # add_surface_elevation(plotaxes,1)
    add_land(plotaxes)
    
    # Bottom surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Internal Surface'
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    # add_surface_elevation(plotaxes,2,bounds=[-300-0.5,-300+0.5])
    add_surface_elevation(plotaxes,2,bounds=internal_surface_limits)
    # add_surface_elevation(plotaxes,2)
    add_land(plotaxes)
    
    # ========================================================================
    #  Depths
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Depths', figno=42)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize':(14,4)}
    
    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Top Layer Depth'
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    add_layer_depth(plotaxes,1)
    add_land(plotaxes)
    
    # Bottom surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Bottom Layer Depth'
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    add_layer_depth(plotaxes,2)
    add_land(plotaxes)
    
    # ========================================================================
    #  Water Speed
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='speed', figno=1)
    plotfigure.show = True
    plotfigure.kwargs = {'figsize':(14,4)}

    # Top layer speed
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents - Top Layer'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_speed(plotaxes,1,bounds=[0.00,0.2])
    add_speed(plotaxes,1,bounds=top_speed_limits)
    # add_speed(plotaxes,1)
    add_land(plotaxes)
    
    # Bottom layer speed
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents - Bottom Layer'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_speed(plotaxes,2,bounds=[0.0,1e-10])
    add_speed(plotaxes,2,bounds=internal_speed_limits)
    # add_speed(plotaxes,2)
    add_land(plotaxes)
    
    # Individual components
    plotfigure = plotdata.new_plotfigure(name='speed_components',figno=401)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(14,14)}
    
    # Top layer
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "X-Velocity - Top Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,1)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_x_velocity(plotaxes,1,bounds=[-1e-10,1e-10])
    add_x_velocity(plotaxes,1)
    add_land(plotaxes)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Y-Velocity - Top Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,2)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_y_velocity(plotaxes,1,bounds=[-0.000125,0.000125])
    add_y_velocity(plotaxes,1)
    add_land(plotaxes)
    
    # Bottom layer
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "X-Velocity - Bottom Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,3)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_x_velocity(plotaxes,2,bounds=[-1e-10,1e-10])
    add_x_velocity(plotaxes,2)
    add_land(plotaxes)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Y-Velocity - Bottom Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,4)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_y_velocity(plotaxes,2,bounds=[-0.8e-6,.8e-6])
    add_y_velocity(plotaxes,2)
    add_land(plotaxes)

    # ========================================================================
    #  Profile Plots
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='profile', figno=4)
    plotfigure.show = False
    
    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = [-2000,20]
    plotaxes.title = "Profile of depth"
    plotaxes.afteraxes = profile_afteraxes
    
    slice_index = 30
    
    # Internal surface
    def bathy_profile(current_data):
        return current_data.x[:,slice_index], b(current_data)[:,slice_index]
    
    def lower_surface(current_data):
        if multilayer_data.init_type == 2:
            return current_data.x[:,slice_index], eta2(current_data)[:,slice_index]
        elif multilayer_data.init_type == 6:
            return current_data.y[slice_index,:], eta2(current_data)[slice_index,:]
        
    
    def upper_surface(current_data):
        if multilayer_data.init_type == 2:
            return current_data.x[:,slice_index], eta1(current_data)[:,slice_index]
        elif multilayer_data.init_type == 6:
            return current_data.y[slice_index,:], eta1(current_data)[slice_index,:]
        
        
    def top_speed(current_data):
        if multilayer_data.init_type == 2:
            return current_data.x[:,slice_index], water_u1(current_data)[:,slice_index]
        elif multilayer_data.init_type == 6:
            return current_data.y[slice_index,:], water_u1(current_data)[slice_index,:]
        
    def bottom_speed(current_data):
        if multilayer_data.init_type == 2:
            return current_data.x[:,slice_index], water_u2(current_data)[:,slice_index]
        elif multilayer_data.init_type == 6:
            return current_data.y[slice_index,:], water_u2(current_data)[slice_index,:]
    
    # Bathy
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = bathy_profile
    plotitem.plot_var = 0
    plotitem.amr_plotstyle = ['-','+','x']
    plotitem.color = 'k'
    plotitem.show = True
    
    # Internal Interface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = lower_surface
    plotitem.plot_var = 7
    plotitem.amr_plotstyle = ['-','+','x']
    plotitem.color = 'b'
    plotitem.show = True
    
    # Upper Interface
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = upper_surface
    plotitem.plot_var = 6
    plotitem.amr_plotstyle = ['-','+','x']
    plotitem.color = (0.2,0.8,1.0)
    plotitem.show = True
    
    # ========================================================================
    #  Combined Profile Plot
    # ========================================================================
    # add_combined_profile_plot(plotdata,0.25,direction='y',figno=120)
    # add_combined_profile_plot(plotdata,0.8,direction='y',figno=121)
    # 
    # add_velocities_profile_plot(plotdata,0.25,direction='y',figno=130)
    # add_velocities_profile_plot(plotdata,0.8,direction='y',figno=131)
    
        
    # ========================================================================
    #  Bathy Profile
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='bathy_profile',figno=20)
    plotfigure.show = False
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.title = "Bathymetry Profile"
    plotaxes.scaled = 'equal'
    
    plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    plotitem.plot_var = b
    plotitem.imshow_cmin = -4000
    plotitem.imshow_cmax = 10
    plotitem.add_colorbar = True
    plotitem.amr_imshow_show = [1,1,1]
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]
    plotitem.show = True
    
    
    # ========================================================================
    # Figure for grids alone
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='grids', figno=11)
    plotfigure.show = False
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = xlimits
    plotaxes.xlimits = ylimits
    plotaxes.title = 'grids'
    plotaxes.afteraxes = pcolor_afteraxes
    plotaxes.scaled = True
    
    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_grid')
    # plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_grid_bgcolor = ['blue','red','green','cyan','yellow']
    plotitem.amr_gridlines_show = [1,1,0,0,0,0]   
    plotitem.amr_gridedges_show = 1
    
    # ========================================================================
    # Figures for momentum
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='x-momentum', figno=13)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'X-Velocity'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    
    # Water
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # # plotitem.plot_var = geoplot.surface
    # plotitem.plot_var = water_u
    # plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    # # plotitem.pcolor_cmin = -1.e-10
    # # plotitem.pcolor_cmax = 1.e-10
    # # plotitem.pcolor_cmin = -2.5 # -3.0
    # # plotitem.pcolor_cmax = 2.5 # 3.0
    # plotitem.add_colorbar = True
    # plotitem.amr_gridlines_show = [0,0,0]
    # plotitem.amr_gridedges_show = [1,1,1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]
    
    plotfigure = plotdata.new_plotfigure(name='y-momentum', figno=14)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Y-Velocity'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    
    # Water
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # # plotitem.plot_var = geoplot.surface
    # plotitem.plot_var = water_v
    # plotitem.pcolor_cmap = colormaps.make_colormap({1.0:'r',0.5:'w',0.0:'b'})
    # # plotitem.pcolor_cmin = -1.e-10
    # # plotitem.pcolor_cmax = 1.e-10
    # # plotitem.pcolor_cmin = -2.5 # -3.0
    # # plotitem.pcolor_cmax = 2.5 # 3.0
    # plotitem.add_colorbar = True
    # plotitem.amr_gridlines_show = [0,0,0]
    # plotitem.amr_gridedges_show = [1,1,1]

    # Land
    add_land(plotaxes)
    
    # ========================================================================
    #  Contour plot for surface
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='contour_surface',figno=15)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(14,4)}
    
    # Set up for axes in this figure:
    
    # Top Surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Top Surface'
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = contour_afteraxes
    add_surface_elevation(plotaxes,plot_type='contour',surface=1,bounds=[-2.5,-1.5,-0.5,0.5,1.5,2.5])
    add_land(plotaxes,plot_type='contour')
    
    # Internal Surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Internal Surface'
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = contour_afteraxes
    add_surface_elevation(plotaxes,plot_type='contour',surface=2,bounds=[-2.5,-1.5,-0.5,0.5,1.5,2.5])
    add_land(plotaxes,plot_type='contour')
    
    # ========================================================================
    #  Contour plot for speed
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='contour_speed',figno=16)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(14,4)}
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Current'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = contour_afteraxes
    
    # Surface
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = water_speed_depth_ave
    plotitem.kwargs = {'linewidths':1}
    # plotitem.contour_levels = [1.0,2.0,3.0,4.0,5.0,6.0]
    plotitem.contour_levels = [0.5,1.5,3,4.5,6.0]
    plotitem.amr_contour_show = [1,1,1]
    plotitem.amr_gridlines_show = [0,0,0]
    plotitem.amr_gridedges_show = [1,1,1]
    plotitem.amr_contour_colors = 'k'
    # plotitem.amr_contour_colors = ['r','k','b']  # color on each level
    # plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.show = True 
    
    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.land
    plotitem.contour_nlevels = 40
    plotitem.contour_min = 0.0
    plotitem.contour_max = 100.0
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.amr_grid_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.gridlines_show = 0
    plotitem.gridedges_show = 0
    plotitem.show = True
    
    # ========================================================================
    #  Vorticity Plot
    # ========================================================================
    # plotfigure = plotdata.new_plotfigure(name='vorticity',figno=17)
    # plotfigure.show = False
    # plotaxes = plotfigure.new_plotaxes()
    # plotaxes.title = "Vorticity"
    # plotaxes.scaled = True
    # plotaxes.xlimits = xlimits
    # plotaxes.ylimits = ylimits
    # plotaxes.afteraxes = pcolor_afteraxes
    # 
    # # Vorticity
    # plotitem = plotaxes.new_plotitem(plot_type='2d_imshow')
    # plotitem.plot_var = 9
    # plotitem.imshow_cmap = plt.get_cmap('PRGn')
    # # plotitem.pcolor_cmap = plt.get_cmap('PuBu')
    # # plotitem.pcolor_cmin = 0.0
    # # plotitem.pcolor_cmax = 6.0
    # plotitem.imshow_cmin = -1.e-2
    # plotitem.imshow_cmax = 1.e-2
    # plotitem.add_colorbar = True
    # plotitem.amr_gridlines_show = [0,0,0]
    # plotitem.amr_gridedges_show = [1]
    # 
    # # Land
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.land
    # plotitem.pcolor_cmap = geoplot.land_colors
    # plotitem.pcolor_cmin = 0.0
    # plotitem.pcolor_cmax = 80.0
    # plotitem.add_colorbar = False
    # plotitem.amr_gridlines_show = [0,0,0]
    
    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0,40.0*3600.0]
    plotaxes.ylimits = top_surface_limits
    plotaxes.title = 'Top Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 6
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    # plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    # plotitem.plot_var = geoplot.gaugetopo
    # plotitem.plotstyle = 'g+'
    
    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                    # create html files of plots?
    plotdata.latex = False                   # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    
