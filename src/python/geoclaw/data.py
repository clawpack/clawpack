#!/usr/bin/env python

"""Classes representing parameters for GeoClaw runs"""

import os

import clawpack.clawutil.data

# Radius of earth in meters.  
# For consistency, should always use this value when needed, e.g.
# in setrun.py or topotools:
Rearth = 6367.5e3  # average of polar and equatorial radii

class GeoClawData(clawpack.clawutil.data.ClawData):
    r"""
    Object containing the basic .

    Note that this data object will write out multiple files.
    """
    def __init__(self):
        super(GeoClawData,self).__init__()

        # GeoClaw physics parameters
        self.add_attribute('gravity',9.8)
        self.add_attribute('earth_radius',Rearth)
        self.add_attribute('coordinate_system',1)
        self.add_attribute('coriolis_forcing',True)
        self.add_attribute('theta_0',45.0)
        self.add_attribute('friction_forcing',True)
        self.add_attribute('manning_coefficient',[0.025])
        self.add_attribute('manning_break',[])

        # GeoClaw algorithm parameters
        self.add_attribute('dry_tolerance',1e-3)
        self.add_attribute('friction_depth',1.0e6)
        self.add_attribute('sea_level',0.0)

        # Refinement behavior
        self.add_attribute('variable_dt_refinement_ratios', False)



    def write(self,data_source='setrun.py'):

        self.open_data_file('geoclaw.data',data_source)

        self.data_write('gravity')
        self.data_write('earth_radius')
        self.data_write('coordinate_system')
        self.data_write('sea_level')
        self.data_write()

        # Forcing terms
        self.data_write('coriolis_forcing')
        if self.coordinate_system == 1 and self.coriolis_forcing:
            self.data_write('theta_0')
        self.data_write('friction_forcing')
        if self.friction_forcing:
            if type(self.manning_coefficient) in [int,float]:
                self.manning_coefficient = [self.manning_coefficient]
            num_manning = len(self.manning_coefficient)
            if len(self.manning_break) != num_manning - 1:
                raise IOError("***manning_break array has wrong length")
            self.data_write(value=num_manning,alt_name='num_manning')
            self.data_write('manning_coefficient')
            self.data_write('manning_break')
            self.data_write('friction_depth')

        self.data_write()

        self.data_write('dry_tolerance',1e-3)
 
        self.close_data_file()



class RefinementData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(RefinementData,self).__init__()

        # Refinement controls
        self.add_attribute('dry_tolerance',1.0e-3)
        self.add_attribute('wave_tolerance',1.0e-1)
        self.add_attribute('speed_tolerance',[1.0e12]*6)
        self.add_attribute('deep_depth',1.0e2)
        self.add_attribute('max_level_deep',3)


    def write(self, data_source="setrun.py"):
        
        # Refinement controls
        self.open_data_file('refinement.data',data_source)
        self.data_write('wave_tolerance')
        if not isinstance(self.speed_tolerance,list):
            self.speed_tolerance = [self.speed_tolerance]
        self.data_write('speed_tolerance')
        self.data_write('deep_depth')
        self.data_write('max_level_deep')
        self.data_write()
        self.data_write('variable_dt_refinement_ratios')

        self.close_data_file()


class TopographyData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(TopographyData,self).__init__()

        # Topography data
        self.add_attribute('test_topography',0)
        self.add_attribute('topofiles',[])
        
        # Jump discontinuity
        self.add_attribute('topo_location',-50e3)
        self.add_attribute('topo_left',-4000.0)
        self.add_attribute('topo_right',-200.0)
        self.add_attribute('topo_angle',0.0)
        
        # Simple oceanic shelf
        self.add_attribute('x0',350e3)
        self.add_attribute('x1',450e3)
        self.add_attribute('x2',480e3)
        self.add_attribute('basin_depth',-3000.0)
        self.add_attribute('shelf_depth',-100.0)
        self.add_attribute('beach_slope',0.008)


    def write(self,data_source='setrun.py'): 

        self.open_data_file('topo.data',data_source)
        self.data_write(name='test_topography',description='(Type topography specification)')
        if self.test_topography == 0:
            ntopofiles = len(self.topofiles)
            self.data_write(value=ntopofiles,alt_name='ntopofiles')
            for tfile in self.topofiles:
                fname = os.path.abspath(tfile[-1])
                self._out_file.write("\n'%s' \n " % fname)
                self._out_file.write("%3i %3i %3i %20.10e %20.10e \n" % tuple(tfile[:-1]))
        elif self.test_topography == 1:
            self.data_write(name='topo_location',description='(Bathymetry jump location)')
            self.data_write(name='topo_left',description='(Depth to left of bathy_location)')
            self.data_write(name='topo_right',description='(Depth to right of bathy_location)')
        elif self.test_topography == 2 or self.test_topography == 3: 
            self.data_write(name='x0',description='(Location of basin end)')
            self.data_write(name='x1',description='(Location of shelf slope end)')
            self.data_write(name='x2',description='(Location of beach slope)')
            self.data_write(name='basin_depth',description='(Depth of basin)')
            self.data_write(name='shelf_depth',description='(Depth of shelf)')
            self.data_write(name='beach_slope',description='(Slope of beach)')
        else:
            raise NotImplementedError("Test topography type %s has not been"
                                        " implemented." % self.test_topography)

        self.close_data_file()


class RefinementData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(RefinementData,self).__init__()

        # Refinement controls
        self.add_attribute('wave_tolerance',1.0e-1)
        self.add_attribute('speed_tolerance',[1.0e12]*6)
        self.add_attribute('deep_depth',1.0e2)
        self.add_attribute('max_level_deep',3)
        self.add_attribute('variable_dt_refinement_ratios',False)


    def write(self,data_source='setrun.py'):
        # Refinement controls
        self.open_data_file('refinement.data',data_source)
        self.data_write('wave_tolerance')
        if not isinstance(self.speed_tolerance,list):
            self.speed_tolerance = [self.speed_tolerance]
        self.data_write('speed_tolerance')
        self.data_write('deep_depth')
        self.data_write('max_level_deep')
        self.data_write()
        self.data_write('variable_dt_refinement_ratios',
                        description="(Set dt refinement ratios automatically)")
        self.close_data_file()



class FixedGridData(clawpack.clawutil.data.ClawData):

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

class FGmaxData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(FGmaxData,self).__init__()
        
        # File name for fgmax points and parameters:
        self.add_attribute('fgmax_files',[])
        self.add_attribute('num_fgmax_val',1)


    def write(self,data_source='setrun.py'):
        self.open_data_file('fgmax.data',data_source)
        num_fgmax_val = self.num_fgmax_val
        if num_fgmax_val not in [1,2,5]:
            raise NotImplementedError( \
                    "Expecting num_fgmax_val in [1,2,5], got %s" % num_fgmax_val)
        self.data_write(value=num_fgmax_val,alt_name='num_fgmax_val')
        num_fgmax_grids = len(self.fgmax_files)
        self.data_write(value=num_fgmax_grids,alt_name='num_fgmax_grids')
        self.data_write()
        for fgmax_file in self.fgmax_files:
            fname = os.path.abspath(fgmax_file)
            self._out_file.write("\n'%s' \n" % fname)
        self.close_data_file()



class DTopoData(clawpack.clawutil.data.ClawData):

    def __init__(self):

        super(DTopoData,self).__init__()
        
        # Moving topograhpy
        self.add_attribute('dtopofiles',[])
        self.add_attribute('dt_max_dtopo', 1.e99)

    def write(self,data_source='setrun.py'):

        # Moving topography settings
        self.open_data_file('dtopo.data',data_source)
        mdtopofiles = len(self.dtopofiles)
        self.data_write(value=mdtopofiles,alt_name='mdtopofiles')
        self.data_write()
        for tfile in self.dtopofiles:
            fname = os.path.abspath(tfile[-1])
            self._out_file.write("\n'%s' \n" % fname)
            self._out_file.write("%3i %3i %3i\n" % tuple(tfile[:-1]))
        self.data_write()
        self.data_write(value=self.dt_max_dtopo,alt_name='dt_max_dtopo')
        self.close_data_file()



class QinitData(clawpack.clawutil.data.ClawData):

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
        if self.qinit_type == 0:
            pass
        elif 0 < self.qinit_type < 5:
            # Check to see if each qinit file is present and then write the data
            for tfile in self.qinitfiles:
                fname = "'%s'" % os.path.abspath(tfile[-1])
                self._out_file.write("\n%s  \n" % fname)
                self._out_file.write("%3i %3i \n" % tuple(tfile[:-1]))
        else:
            raise ValueError("Invalid qinit_type parameter %s." % self.qinit_type)
        self.close_data_file()
