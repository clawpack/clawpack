
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 
    
def setplot(plotdata,  bathy_location=0.15,  bathy_angle=0.0,  
                       bathy_left=-1.0,      bathy_right=-0.2):
    """Setup the plotting data objects.

    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    returns plotdata object

    """ 

    import os

    import numpy as np
    import matplotlib.pyplot as plt

    from clawpack.visclaw import geoplot, gaugetools

    import clawpack.clawutil.data as clawutil
    import clawpack.amrclaw.data as amrclaw
    import clawpack.geoclaw.data

    import clawpack.geoclaw.multilayer.plot as ml_plot

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir,'claw.data'))
    amrdata = amrclaw.AmrclawInputData(clawdata)
    amrdata.read(os.path.join(plotdata.outdir,'amr.data'))
    geodata = clawpack.geoclaw.data.GeoClawData()
    geodata.read(os.path.join(plotdata.outdir,'geoclaw.data'))
    multilayer_data = clawpack.geoclaw.data.MultilayerData()
    multilayer_data.read(os.path.join(plotdata.outdir,'multilayer.data'))

    def transform_c2p(x,y,x0,y0,theta):
        return ((x+x0)*np.cos(theta) - (y+y0)*np.sin(theta),
                (x+x0)*np.sin(theta) + (y+y0)*np.cos(theta))

    def transform_p2c(x,y,x0,y0,theta):
        return ( x*np.cos(theta) + y*np.sin(theta) - x0,
                -x*np.sin(theta) + y*np.cos(theta) - y0)
        
    
    # Setup bathymetry reference lines
    with open(os.path.join(plotdata.outdir,"bathy_geometry.data"), 'r') as bathy_geometry_file:
        bathy_location = float(bathy_geometry_file.readline())
        bathy_angle = float(bathy_geometry_file.readline())
    x = [0.0,0.0]
    y = [0.0,1.0]
    x1,y1 = transform_c2p(x[0],y[0],bathy_location,
                                0.0,bathy_angle)
    x2,y2 = transform_c2p(x[1],y[1],bathy_location,
                                0.0,bathy_angle)
    
    if abs(x1 - x2) < 10**-3:
        x = [x1, x1]
        y = [clawdata.lower[1], clawdata.upper[1]]
    else:
        m = (y1 - y2) / (x1 - x2)
        x[0] = (clawdata.lower[1] - y1) / m + x1
        y[0] = clawdata.lower[1]
        x[1] = (clawdata.upper[1] - y1) / m + x1
        y[1] = clawdata.upper[1]
    ref_lines = [((x[0],y[0]),(x[1],y[1]))]

    plotdata.clearfigures()
    
    plotdata.save_frames = False
    
    # ========================================================================
    #  Generic helper functions
    # ========================================================================
    def pcolor_afteraxes(current_data):
        bathy_ref_lines(current_data)
        gauge_locations(current_data)
        
    def contour_afteraxes(current_data):
        # gauge_locations(current_data)
        # m_to_km_labels()
        plt.hold(True)
        pos = -80.0 * (23e3 / 180) + 500e3 - 5e3
        plt.plot([pos,pos],[-300e3,300e3],'b',[pos-5e3,pos-5e3],[-300e3,300e3],'y')
        plt.hold(False)
        wind_contours(current_data)
        bathy_ref_lines(current_data)
        
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
        plt.hold(True)
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos=gaugenos, format_string='kx', add_labels=True)
        plt.hold(False)


    # ========================================================================
    # Axis limits
    #xlimits = [amrdata.xlower,amrdata.xupper]
    xlimits = [-0.5,0.5]
    #ylimits = [amrdata.ylower,amrdata.yupper]
    ylimits = [-0.5,0.5]
    eta = [multilayer_data.eta[0],multilayer_data.eta[1]]
    top_surface_limits = [eta[0]-0.03,eta[0]+0.03]
    internal_surface_limits = [eta[1]-0.015,eta[1]+0.015]
    # top_surface_limits = [eta[0]-0.3,eta[0]+0.3]
    # internal_surface_limits = [eta[1]-0.15,eta[1]+0.15]
    top_speed_limits = [0.0,0.1]
    internal_speed_limits = [0.0,0.03]
    
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
    ml_plot.add_surface_elevation(plotaxes,1,bounds=top_surface_limits)
    # ml_plot.add_surface_elevation(plotaxes,1,bounds=[-0.06,0.06])
    # ml_plot.add_surface_elevation(plotaxes,1)
    ml_plot.add_land(plotaxes)
    
    # Bottom surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Internal Surface'
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    # ml_plot.add_surface_elevation(plotaxes,2,bounds=[-300-0.5,-300+0.5])
    ml_plot.add_surface_elevation(plotaxes,2,bounds=internal_surface_limits)
    # ml_plot.add_surface_elevation(plotaxes,2)
    ml_plot.add_land(plotaxes)
    
    # ========================================================================
    #  Depths
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Depths', figno=42)
    plotfigure.show = False
    plotfigure.kwargs = {'figsize':(14,4)}
    
    # Top surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Top Layer Depth'
    plotaxes.axescmd = 'subplot(1,2,1)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    ml_plot.add_layer_depth(plotaxes,1,bounds=[-0.1,1.1])
    ml_plot.add_land(plotaxes)
    
    # Bottom surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Bottom Layer Depth'
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = pcolor_afteraxes
    ml_plot.add_layer_depth(plotaxes,2,bounds=[-0.1,0.7])
    ml_plot.add_land(plotaxes)
    
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
    ml_plot.add_speed(plotaxes,1,bounds=top_speed_limits)
    # add_speed(plotaxes,1)
    ml_plot.add_land(plotaxes)
    
    # Bottom layer speed
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Currents - Bottom Layer'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_speed(plotaxes,2,bounds=[0.0,1e-10])
    ml_plot.add_speed(plotaxes,2,bounds=internal_speed_limits)
    # add_speed(plotaxes,2)
    ml_plot.add_land(plotaxes)
    
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
    ml_plot.add_x_velocity(plotaxes,1)
    ml_plot.add_land(plotaxes)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Y-Velocity - Top Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,2)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_y_velocity(plotaxes,1,bounds=[-0.000125,0.000125])
    ml_plot.add_y_velocity(plotaxes,1)
    ml_plot.add_land(plotaxes)
    
    # Bottom layer
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "X-Velocity - Bottom Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,3)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_x_velocity(plotaxes,2,bounds=[-1e-10,1e-10])
    ml_plot.add_x_velocity(plotaxes,2)
    ml_plot.add_land(plotaxes)
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Y-Velocity - Bottom Layer"
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.axescmd = 'subplot(2,2,4)'
    plotaxes.afteraxes = pcolor_afteraxes
    # add_y_velocity(plotaxes,2,bounds=[-0.8e-6,.8e-6])
    ml_plot.add_y_velocity(plotaxes,2)
    ml_plot.add_land(plotaxes)

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
    plotitem.plot_var = ml_plot.b
    plotitem.imshow_cmin = -4000
    plotitem.imshow_cmax = 10
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    plotitem.show = True
    
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
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.amr_patchedges_show = [1,1,1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.show = True
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 80.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
    
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
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.amr_patchedges_show = [1,1,1]

    # Land
    ml_plot.add_land(plotaxes)
    
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
    ml_plot.add_surface_elevation(plotaxes,plot_type='contour',surface=1,bounds=[-2.5,-1.5,-0.5,0.5,1.5,2.5])
    ml_plot.add_land(plotaxes,plot_type='contour')
    
    # Internal Surface
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = 'Internal Surface'
    plotaxes.axescmd = 'subplot(1,2,2)'
    plotaxes.scaled = True
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = contour_afteraxes
    ml_plot.add_surface_elevation(plotaxes,plot_type='contour',surface=2,bounds=[-2.5,-1.5,-0.5,0.5,1.5,2.5])
    ml_plot.add_land(plotaxes,plot_type='contour')
    
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
    plotitem.plot_var = ml_plot.water_speed_depth_ave
    plotitem.kwargs = {'linewidths':1}
    # plotitem.contour_levels = [1.0,2.0,3.0,4.0,5.0,6.0]
    plotitem.contour_levels = [0.5,1.5,3,4.5,6.0]
    plotitem.amr_contour_show = [1,1,1]
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.amr_patchedges_show = [1,1,1]
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
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = 0
    plotitem.amr_patchedges_show = 0
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
    # plotitem.amr_celledges_show = [0,0,0]
    # plotitem.amr_patchedges_show = [1]
    # 
    # # Land
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.land
    # plotitem.pcolor_cmap = geoplot.land_colors
    # plotitem.pcolor_cmin = 0.0
    # plotitem.pcolor_cmax = 80.0
    # plotitem.add_colorbar = False
    # plotitem.amr_celledges_show = [0,0,0]
    
    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    def plot_regression_gauges(cd, plot_var=(3, 7)):
        import clawpack.pyclaw.gauges as gauges
        import numpy

        fig = plt.gcf()
        gauge = gauges.GaugeSolution(cd.gaugeno, path="./regression_data")
        for (i, n) in enumerate(plot_var):
            fig.axes[i].plot(gauge.t, gauge.q[n, :], 'k-')

    # Top
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0,1.0]
    # plotaxes.ylimits = top_surface_limits
    plotaxes.title = 'Top Surface'
    plotaxes.axescmd = "subplot(1, 2, 1)"

    # Plot surfaces:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'ro'

    # Bottom
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0,1.0]
    plotaxes.title = 'Bottom Surface'
    plotaxes.axescmd = "subplot(1, 2, 2)"
    plotaxes.afteraxes = plot_regression_gauges

    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 7
    plotitem.plotstyle = 'ro'

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

    
