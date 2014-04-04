# encoding: utf-8
""" 
Module to setup storm surge parameters
""" 

import clawpack.clawutil.data as clawdata

# Simple hurricane data format
class SurgeData(clawdata.ClawData):

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
        
        # self.data_write('wind_tolerance',description='(Wind speed tolerance)')
        # self.data_write('pressure_tolerance',description="(Pressure source term tolerance)")
        # self.data_write()
                
        self.data_write('wind_refine',description='(Refinement ratios)')
        self.data_write('R_refine',description="(Refinement ratios)")
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
            self.close_data_file()
            raise ValueError("Invalid storm type %s." % self.storm_type)

        self.close_data_file()



class FrictionData(clawdata.ClawData):
    r"""Data class representing variable friction parameters and data sources"""

    def __init__(self):

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