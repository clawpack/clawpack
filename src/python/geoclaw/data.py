#!/usr/bin/env python

"""Classes representing parameters for GeoClaw runs"""

import clawpack.clawutil.clawdata

class GeoclawInputData(clawpack.clawutil.clawdata.ClawData):
    r"""
    Object that will be written out to the various GeoClaw data files.

    Note that this data object will write out multiple files.
    """
    def __init__(self, num_dim):
        super(GeoclawInputData,self).__init__()

        # GeoClaw physics parameters
        self.add_attribute('gravity',9.8)
        self.add_attribute('rho',1025.0)
        self.add_attribute('earth_radius',6367500.0)
        self.add_attribute('coordinate_system',1)
        self.add_attribute('coriolis_forcing',True)
        self.add_attribute('theta_0',45.0)
        self.add_attribute('friction_forcing',True)
        self.add_attribute('manning_coefficient',0.025)
        self.add_attribute('wet_manning_coefficient',None)
        self.add_attribute('dry_manning_coefficient',None)

        # GeoClaw algorithm parameters
        self.add_attribute('friction_depth',1.0e6)
        self.add_attribute('sea_level',0.0)
        self.add_attribute('variable_dt_refinement_ratios',False)


    def write(self,data_source='setrun.py'):

        self.open_data_file('physics.data',data_source)

        self.data_write('gravity')
        self.data_write('rho',description='(Density of water)')
        self.data_write('earth_radius')
        self.data_write('coordinate_system',description="(Coorindat sytem to use, lat-long etc.)")
        self.data_write('sea_level',description='(Nominal sea level)')
        self.data_write()
        self.data_write('coriolis_forcing')

        if self.coordinate_system == 1 and self.coriolis_forcing:
            self.data_write('theta_0')
        self.data_write('friction_forcing')    
        # Write out wet/dry coefficient if provided, otherwise write out the 
        # generic version    
        if self.wet_manning_coefficient is None:
            self.wet_manning_coefficient = self.manning_coefficient
        self.data_write('wet_manning_coefficient',description='(Manning coefficient used when initially wet)')    
        if self.dry_manning_coefficient is None:
            self.dry_manning_coefficient = self.manning_coefficient    
        self.data_write('dry_manning_coefficient',description='(Manning coefficient used when initially dry)')
        self.data_write('friction_depth')
        self.data_write()
        
        self.data_write('dry_tolerance')
        self.data_write('variable_dt_refinement_ratios')


class TopographyData(clawpack.clawutil.clawdata.ClawData):

    def __init__(self):

        super(TopographyData,self).__init__()
        
        # Topography data
        self.add_attribute('topo_type',0)
        self.add_attribute('topofiles',[])
        
        self.add_attribute('topo_location',-50e3)
        self.add_attribute('topo_left',-4000.0)
        self.add_attribute('topo_right',-200.0)
        self.add_attribute('topo_angle',0.0)
        
        self.add_attribute('x0',350e3)
        self.add_attribute('x1',450e3)
        self.add_attribute('x2',480e3)
        self.add_attribute('basin_depth',-3000.0)
        self.add_attribute('shelf_depth',-100.0)
        self.add_attribute('beach_slope',0.008)


    def write(self,data_source='setrun.py'):        

        # Topography data
        self.open_data_file('topo.data',data_source)
        self.data_write(name='topo_type',description='(Type topography specification)')
        if self.topo_type == 0:
            ntopofiles = len(self.topofiles)
            self.data_write(value=ntopofiles,alt_name='ntopofiles')
            for tfile in self.topofiles:
                try:
                    fname = os.path.abspath(tfile[-1])
                except:
                    print "*** Error: file not found: ",tfile[-1]
                    raise MissingFile("file not found")
                self._out_file.write("\n'%s' \n " % fname)
                self._out_file.write("%3i %3i %3i %20.10e %20.10e \n" % tuple(tfile[:-1]))
        elif self.topo_type == 1:
            self.data_write(name='topo_location',description='(Bathymetry jump location)')
            self.data_write(name='topo_left',description='(Depth to left of bathy_location)')
            self.data_write(name='topo_right',description='(Depth to right of bathy_location)')
        elif self.topo_type == 2 or self.topo_type == 3: 
            self.data_write(name='x0',description='(Location of basin end)')
            self.data_write(name='x1',description='(Location of shelf slope end)')
            self.data_write(name='x2',description='(Location of beach slope)')
            self.data_write(name='basin_depth',description='(Depth of basin)')
            self.data_write(name='shelf_depth',description='(Depth of shelf)')
            self.data_write(name='beach_slope',description='(Slope of beach)')
        else:
            raise NotImplemented("Topography type %s has not been implemented." % topo_type)    
        self.close_data_file()


class RefinementData(clawpack.clawutil.clawdata.ClawData):

    def __init__(self):

        super(RefinementData,self).__init__()

        # Refinement controls
        self.add_attribute('dry_tolerance',1.0e-3)
        self.add_attribute('wave_tolerance',1.0e-1)
        self.add_attribute('speed_tolerance',[1.0e12]*6)
        self.add_attribute('deep_depth',1.0e2)
        self.add_attribute('max_level_deep',3)


    def write(self,data_source='setrun.py'):
        # Refinement controls
        self.open_data_file('refinement.data',data_source)
        self.data_write('wave_tolerance')
        if not isinstance(self.speed_tolerance,list):
            self.speed_tolerance = [self.speed_tolerance]
            print "Warning: Need len(speed_tolerance) >= mxnest"
        self.data_write('speed_tolerance')
        self.data_write('deep_depth')
        self.data_write('max_level_deep')
        self.data_write()
        

class FixedGridData(clawpack.clawutil.clawdata.ClawData):

    def __init__(self):

        super(FixedGridData,self).__init__()
        
        # Fixed Grids
        self.add_attribute('fixedgrids',[])


    def write(self,data_source='setrun.py'):
        # Fixed grid settings
        self.open_data_file('fixed_grids.data',data_source)
        nfixedgrids = len(self.fixedgrids)
        self.data_write(value=nfixedgrids,alt_name='nfixedgrids')
        self.data_write()
        for fixedgrid in self.fixedgrids:
            self._out_file.write(11*"%g  " % tuple(fixedgrid) +"\n")
        self.close_data_file()


class DTopoData(clawpack.clawutil.clawdata.ClawData):

    def __init__(self):

        super(DTopoData,self).__init__()
        
        # Moving topograhpy
        self.add_attribute('dtopofiles',[])

    def write(self,data_source='setrun.py'):
        # Moving topography settings
        self.open_data_file('dtopo.data',data_source)
        mdtopofiles = len(self.dtopofiles)
        self.data_write(value=mdtopofiles,alt_name='mdtopofiles')
        self.data_write()
        for tfile in self.dtopofiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                # print "*** Error: file not found: ",tfile[-1]
                raise MissingFile("file not found")
            self._out_file.write("\n%s \n" % fname)
            self._out_file.write("%3i %3i %3i\n" % tuple(tfile[:-1]))
        self.close_data_file()


class QinitData(clawpack.clawutil.clawdata.ClawData):

    def __init__(self):

        super(QinitData,self).__init__()
        
        # Qinit data
        self.add_attribute('qinit_type',0)
        self.add_attribute('qinitfiles',[])    

    def write(self,data_source='setrun.py'):
        # Initial perturbation
        self.open_data_file('qinit.data',data_source)
        self.data_write('qinit_type')

        # Perturbation requested
        self.data_write()
        # Check to see if each qinit file is present and then write the data
        for tfile in self.qinitfiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                raise MissingFile("file not found")
            self._out_file.write("\n%s  \n" % fname)
            self._out_file.write("%3i %3i \n" % tuple(tfile[:-1]))
        self.close_data_file()