#!/usr/bin/env python

"""

Classes representing parameters for GeoClaw runs

:Classes:

 - GeoClawData
 - RefinementData
 - TopographyData
 - FixedGridData
 - FGmaxData
 - DTopoData
 - QinitData

:Constants:

 - Rearth - Radius of earth in meters
 - DEG2RAD factor to convert degrees to radians
 - RAD2DEG factor to convert radians to degrees
 - LAT2METER factor to convert degrees in latitude to meters
"""

import os
import numpy
import clawpack.clawutil.data

# Radius of earth in meters.  
# For consistency, should always use this value when needed, e.g.
# in setrun.py or topotools:
Rearth = 6367.5e3  # average of polar and equatorial radii

DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi
LAT2METER = Rearth * DEG2RAD

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

        self.data_write('dry_tolerance')
 
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


    def read(self, path, force=False):
        r"""Read a dtopography data file."""

        print self.dtopofiles

        with open(os.path.abspath(path), 'r') as data_file:

            file_name = None
            
            # Forward to first parameter
            for line in data_file:

                # Regular parameter setting
                if "=:" in line:
                    value, tail = line.split("=:")
                    varname = tail.split()[0]

                    if varname == "mdtopofiles":
                        num_dtopo_files = int(value)
                    elif varname == "dt_max_dtopo":
                        self.dt_max_dtopo = float(value)

                # Assume this is the second line of a record, complete and add
                # to dtopofiles list
                elif file_name is not None:
                    base_values = [int(value) for value in line.split()]
                    base_values.append(file_name)
                    self.dtopofiles.append(base_values)
                    file_name = None

                # Non-empty line, assume start of dtopo_file record
                elif line[0] == "'":
                    file_name = line.strip()[1:-1]

        # Check to make sure we have all the dtopofiles
        if len(self.dtopofiles) != num_dtopo_files:
            raise IOError("The number of dtopo files specified does not equal ",
                          "the number found.")


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


# Storm data
class SurgeData(clawpack.clawutil.data.ClawData):
    r"""Data object describing storm surge related parameters"""

    def __init__(self):
        super(SurgeData,self).__init__()

        # Physics parameters
        self.add_attribute('rho_air',1.15) # Density of air
        self.add_attribute('ambient_pressure',101.3e3) # Nominal atmos pressure

        # Source term controls
        self.add_attribute('wind_forcing',False)
        self.add_attribute('drag_law',1)
        self.add_attribute('pressure_forcing',False)
        
        # Source term algorithm parameters
        # self.add_attribute('wind_tolerance',1e-6)
        # self.add_attribute('pressure_tolerance',1e-4) # Pressure source term tolerance

        # AMR parameters
        self.add_attribute('wind_refine',[20.0,40.0,60.0])
        self.add_attribute('R_refine',[60.0e3,40e3,20e3])
        
        # Storm parameters
        self.add_attribute("storm_type",0) # Type of storm
        self.add_attribute("landfall",0.0)

        # Storm type 1 - Read in file track
        self.add_attribute("storm_file",'./storm.data')

        # Storm type 2 - Idealized Hurricane, these match hurricane Tracy
        self.add_attribute("ramp_up_t",12*60.0**2) # Ramp up time for hurricane
        self.add_attribute('velocity',(5.0,0.0))  # Speed of hurricane
        self.add_attribute('R_eye_init',(0.0,0.0))     # Initial position
        self.add_attribute('A',23.0)        # Hurricane model fit parameter
        self.add_attribute('B',1.5)
        self.add_attribute('Pc',950.0)      # Pressure in the eye of the hurricane 

        # Storm type 3 - Stommel wind field
        self.add_attribute('stommel_wind',1.0)

        # Algorithm parameters
        self.add_attribute("wind_index", 5)
        self.add_attribute("pressure_index", 7)

        
    def write(self,out_file='./surge.data',data_source="setrun.py"):
        """Write out the data file to the path given"""

        # print "Creating data file %s" % out_file
        self.open_data_file(out_file,data_source)

        self.data_write('rho_air',description="(Density of air)")
        self.data_write('ambient_pressure',description="(Nominal atmospheric pressure)")
        self.data_write()
        
        self.data_write('wind_forcing',description='(Wind source term used)')
        self.data_write('drag_law',description='(Type of drag law to use)')
        self.data_write('pressure_forcing',description="(Pressure source term used)")
        self.data_write()

        self.data_write("wind_index", description="(Index into aux array for wind (size 2))")
        self.data_write("pressure_index", description="(Index into aux array for pressure (size 1))")
        self.data_write()

        if isinstance(self.wind_refine, bool):
            if not self.wind_refine:
                self.data_write('wind_refine', value=False, description='(Refinement ratios)')
        elif isinstance(self.wind_refine, type(None)):
            self.data_write('wind_refine', value=False, description='(Refinement ratios)')
        else:
            self.data_write('wind_refine',description='(Refinement ratios)')
        if isinstance(self.R_refine, bool):
            if not self.R_refine:
                self.data_write('R_refine', value=False, description='(Refinement ratios)')
        elif isinstance(self.R_refine, type(None)):
            self.data_write('R_refine', value=False, description='(Refinement ratios)')
        else:
            self.data_write('R_refine',description='(Refinement ratios)')
        self.data_write()
        
        self.data_write("storm_type",description='(Storm specification type)')
        self.data_write('landfall',description="(Landfall time of storm)")
        self.data_write('storm_file',description="(Location of storm data)")

        if self.storm_type == 0 or self.storm_type == 1:
            pass 
        elif self.storm_type == 2:
            # Open another data file called stored in storm_file and write the 
            # following parameters to it
            self.open_data_file(self.storm_file)
            self.data_write("ramp_up_t",description="(Ramp up time for wind field)")
            self.data_write('velocity',description="(Speed of storm)")
            self.data_write('R_eye_init',description="(Initial position of storm)")
            self.data_write('A',description="(Hurricane model fit parameter)")
            self.data_write('B')
            self.data_write('Pc',description="(Pressure in the eye of the hurricane)")
        elif self.storm_type == 3:
            # Open another data file called stored in storm_file and write the 
            # following parameters to it
            self.open_data_file(self.storm_file)
            self.data_write("stommel_wind",description="(Amplitude of Stommel wind)")
        else:
            raise ValueError("Invalid storm type %s." % self.storm_type)

        self.close_data_file()



class FrictionData(clawpack.clawutil.data.ClawData):
    r"""Data class representing complex variable friction"""

    def __init__(self):
        r""""""
        
        super(FrictionData, self).__init__()

        # Variable friction support
        self.add_attribute('variable_friction',False)

        # Region support
        self.add_attribute('friction_regions',[])

        # File support
        self.add_attribute('friction_files',[])


    def write(self, out_file='friction.data', data_source='setrun.py'):

        self.open_data_file(out_file,data_source)

        self.data_write('variable_friction',description="(method for setting variable friction)")
        self.data_write()
        if self.variable_friction:
            # Region based friction
            self.data_write(value=len(self.friction_regions),
                            alt_name='num_friction_regions',
                            description="(Friction Regions)")
            self.data_write()
            for region in self.friction_regions:
                self.data_write(value=region[0],alt_name="lower")
                self.data_write(value=region[1],alt_name="upper")
                self.data_write(value=region[2],alt_name="depths")
                self.data_write(value=region[3],alt_name="manning_coefficients")
                self.data_write()

            # File based friction
            self.data_write(value=len(self.friction_files),
                            alt_name='num_friction_files')
            for friction_file in self.friction_files:
                self._out_file.write("'%s' %s\n " % friction_file)

        self.close_data_file()



class MultilayerData(clawpack.clawutil.data.ClawData):
    r"""
    Multilayer SWE data object

    """

    def __init__(self):
        super(MultilayerData, self).__init__()

        # Physics parameters
        self.add_attribute('num_layers',2)
        self.add_attribute('rho',[1025.0,1028.0])
        self.add_attribute('eta',[0.0,-200.0])
        
        # Algorithm parameters
        self.add_attribute('eigen_method',4)
        self.add_attribute('inundation_method',2)
        self.add_attribute('check_richardson',True)
        self.add_attribute('richardson_tolerance',0.95)

        # Need to adjust refinement module for this, dry_limit is in geodata
        self.add_attribute('wave_tolerance',[1e-1,2e-1])
        self.add_attribute('dry_limit',False)
    

    def write(self,out_file='./multilayer.data',datasource="setrun.py"):
        
        self.open_data_file(out_file, datasource)
        
        self.data_write('num_layers', description='(Number of layers)')
        self.data_write('rho',description='(Densities of layers)')
        self.data_write('eta',description='(Initial top surface of each layer)')
        self.data_write(None)
        self.data_write('check_richardson',description="(Check Richardson number)")
        self.data_write('richardson_tolerance',description='(Tolerance for Richardson number)')
        self.data_write('eigen_method',description='(Method for calculating eigenspace)')
        self.data_write('inundation_method',description='(Method for calculating inundation eigenspace)')
        
        # self.data_write('wave_tolerance',description='(Tolerance for wave height refinement)')
        # self.data_write('dry_limit',description='(Turn off limiting when near a dry state)')
        
        self.close_data_file()