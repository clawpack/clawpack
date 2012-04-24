#!/usr/bin/env python
# encoding: utf-8
""" 
Module to setup hurrican parameters
""" 

import numpy as np

import clawutil.oldclawdata as data

# Simple hurricane data format
class HurricaneData(data.Data):
    def __init__(self,ramp_up_time):
        super(HurricaneData,self).__init__()

        # Wind source term controls
        self.add_attribute('wind_src',True)
        self.add_attribute('pressure_src',True)
        
        # Source term algorithm parameters
        self.add_attribute('wind_tolerance',1e-6)
        self.add_attribute('pressure_tolerance',1e-4) # Pressure source term tolerance

        # AMR parameters for follwing the storm
        self.add_attribute('max_wind_nest',0) # Wind strength based refinement
        self.add_attribute('wind_refine',[20.0,40.0,60.0])
        self.add_attribute('max_R_nest',3) # Hurricane location based refinement
        self.add_attribute('R_refine',[60.0e3,40e3,20e3])
        
        # Physics parameters
        self.add_attribute('rho_air',1.15) # Density of air
        self.add_attribute('theta_0',30.0)  # Center of beta-plane approximation  
        
        # Wind parameters
        self.add_attribute("wind_type",1) # Type of wind
        
        # Read in hurricane track data
        self.add_attribute("hurricane_track_file",'./hurricane_track.data')
        
        # Idealized Hurricane, these match hurricane Tracy as a default
        self.add_attribute("ramp_up_t",ramp_up_time) # Ramp up time for hurricane
        self.add_attribute('hurricane_velocity',(5.0,0.0))  # Speed of hurricane
        self.add_attribute('R_eye_init',(0.0,0.0))     # Initial position
        self.add_attribute('A',23.0)        # Hurricane model fit parameter
        self.add_attribute('B',1.5)     
        self.add_attribute('Pn',1005.0)     # Nominal atmospheric pressure     
        self.add_attribute('Pc',950.0)      # Pressure in the eye of the hurricane  
        
        # Stommel wind field
        # A is already defined above
        
    def write(self,out_file='./hurricane.data',datasource="setrun.py"):
        """Write out the data file to the path given"""

        print "Creating data file %s" % out_file
        out_file = data.open_datafile(out_file)
        
        data.data_write(out_file,self,'wind_src','(Wind source term used)')
        data.data_write(out_file,self,'pressure_src',"(Pressure source term used)")
        data.data_write(out_file,self,None)
        
        data.data_write(out_file,self,'wind_tolerance','(Wind speed tolerance)')
        data.data_write(out_file,self,'pressure_tolerance',"(Pressure source term tolerance)")
        data.data_write(out_file,self,None)
                
        data.data_write(out_file,self,'max_wind_nest','(Wind strength based refinement)')
        data.data_write(out_file,self,'wind_refine','(Refinement ratios)')
        data.data_write(out_file,self,'max_R_nest',"(Hurricane location based refinement)")
        data.data_write(out_file,self,'R_refine',"(Refinement ratios)")
        data.data_write(out_file,self,None)
        
        data.data_write(out_file,self,'rho_air',"(Density of air)")
        data.data_write(out_file,self,"theta_0","(Center of beta-plane approximation)")
        data.data_write(out_file,self,None)
        
        data.data_write(out_file,self,"wind_type",'(Type of wind field)')
        if self.wind_type == 0:
            data.data_write(out_file,self,"hurricane_track_file",'(Data file containing track)')
        elif self.wind_type == 1:
            data.data_write(out_file,self,"ramp_up_t","(Ramp up time for wind field)")
            data.data_write(out_file,self,'hurricane_velocity',"(Speed of moving wind field)")
            data.data_write(out_file,self,'R_eye_init',"(Initial position of wind field)")
            data.data_write(out_file,self,'A',"(Hurricane model fit parameter)")
            data.data_write(out_file,self,'B')     
            data.data_write(out_file,self,'Pn',"(Nominal atmospheric pressure)")
            data.data_write(out_file,self,'Pc',"(Pressure in the eye of the hurricane)")
        elif self.wind_type == 2:
            data.data_write(out_file,self,'A',"(Wind maximum amplitude)")

        data.data_write(out_file,self,None)
        
        out_file.close()