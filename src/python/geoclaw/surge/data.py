# encoding: utf-8
""" 
Module to setup storm surge parameters
""" 

import numpy as np

import clawpack.clawutil.oldclawdata as data

# Simple hurricane data format
class StormData(data.Data):

    def __init__(self):
        super(StormData,self).__init__()
        
        # Physics parameters
        self.add_attribute('rho_air',1.15) # Density of air
        self.add_attribute('ambient_pressure',101.3e3) # Nominal atmos pressure

        # Source term controls
        self.add_attribute('wind_forcing',True)
        self.add_attribute('pressure_forcing',True)
        
        # Source term algorithm parameters
        self.add_attribute('wind_tolerance',1e-6)
        self.add_attribute('pressure_tolerance',1e-4) # Pressure source term tolerance

        # AMR parameters
        self.add_attribute('wind_refine',[20.0,40.0,60.0])
        self.add_attribute('R_refine',[60.0e3,40e3,20e3])
        
        # Storm parameters
        self.add_attribute("storm_type",1) # Type of storm

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

        
    def write(self,out_file='./surge.data',datasource="setrun.py"):
        """Write out the data file to the path given"""

        # print "Creating data file %s" % out_file
        data_file = data.open_datafile(out_file)

        data.data_write(data_file,self,'rho_air',"(Density of air)")
        data.data_write(data_file,self,'ambient_pressure',"(Nominal atmospheric pressure)")
        data.data_write(data_file,self,None)
        
        data.data_write(data_file,self,'wind_forcing','(Wind source term used)')
        data.data_write(data_file,self,'pressure_forcing',"(Pressure source term used)")
        data.data_write(data_file,self,None)
        
        data.data_write(data_file,self,'wind_tolerance','(Wind speed tolerance)')
        data.data_write(data_file,self,'pressure_tolerance',"(Pressure source term tolerance)")
        data.data_write(data_file,self,None)
                
        data.data_write(data_file,self,'wind_refine','(Refinement ratios)')
        data.data_write(data_file,self,'R_refine',"(Refinement ratios)")
        data.data_write(data_file,self,None)
        
        data.data_write(data_file,self,"storm_type",'(Storm specification type)')
        data.data_write(data_file,self,'storm_file',"(Location of storm data)")
        if self.storm_type == 0 or self.storm_type == 1:
            pass 
        elif self.storm_type == 2:
            # Open another data file called stored in storm_file and write the 
            # following parameters to it
            storm_file = data.open_datafile(self.storm_file)
            data.data_write(storm_file,self,"ramp_up_t","(Ramp up time for wind field)")
            data.data_write(storm_file,self,'velocity',"(Speed of storm)")
            data.data_write(storm_file,self,'R_eye_init',"(Initial position of storm)")
            data.data_write(storm_file,self,'A',"(Hurricane model fit parameter)")
            data.data_write(storm_file,self,'B')
            data.data_write(storm_file,self,'Pc',"(Pressure in the eye of the hurricane)")
        elif self.storm_type == 3:
            # Open another data file called stored in storm_file and write the 
            # following parameters to it
            storm_file = data.open_datafile(self.storm_file)
            data.data_write(storm_file,self,"stommel_wind","(Amplitude of Stommel wind)")
        else:
            raise Exception("Invalid storm type %s." % storm_type)
            

        data.data_write(data_file,self,None)
        
        data_file.close()