#!/usr/bin/env python

"""Classes representing parameters for GeoClaw Multilayer runs"""

import clawpack.clawutil.data

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
        self.add_attribute('richardson_tolerance',0.95)
        self.add_attribute('wave_tolerance',[1e-1,2e-1])
        self.add_attribute('dry_limit',False)
    

    def write(self,out_file='./multilayer.data',datasource="setrun.py"):
        
        self.open_data_file(out_file, datasource)
        
        self.data_write(out_file,self,'layers','(Number of layers)')
        self.data_write(out_file,self,'rho','(Densities of layers)')
        self.data_write(out_file,self,'eta','(Initial top surface of each layer)')
        self.data_write(out_file,self,None)
        self.data_write(out_file,self,'eigen_method','(Method for calculating eigenspace)')
        self.data_write(out_file,self,'inundation_method','(Method for calculating inundation eigenspace)')
        self.data_write(out_file,self,'richardson_tolerance','(Tolerance for Richardson number)')
        self.data_write(out_file,self,'wave_tolerance','(Tolerance for wave height refinement)')
        self.data_write(out_file,self,'dry_limit','(Turn off limiting when near a dry state)')
        
        self.close_data_file()


