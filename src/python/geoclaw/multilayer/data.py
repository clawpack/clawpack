#!/usr/bin/env python

"""Classes representing parameters for GeoClaw Multilayer runs"""

import os

import clawpack.clawutil.data
import clawpack.geoclaw.data

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
        # self.add_attribute('wave_tolerance',[1e-1,2e-1])
        # self.add_attribute('dry_limit',False)
    

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

class QinitMultilayerData(clawpack.geoclaw.data.QinitData):
    r"""
    Modified Qinit data object for multiple layers

    """

    def __init__(self):
        super(QinitMultilayerData, self).__init__()

        # Test qinit data > 4
        self.add_attribute("init_location", [0.0, 0.0])
        self.add_attribute("wave_family", 1)
        self.add_attribute("epsilon", 0.02)
        self.add_attribute("angle", 0.0)
        self.add_attribute("sigma", 0.02)

    def write(self, data_source='setrun.py'):

        # Initial perturbation
        self.open_data_file('qinit.data',data_source)
        self.data_write('qinit_type')

        # Perturbation requested
        if self.qinit_type == 0:
            pass
        elif 0 < self.qinit_type < 5:
            # Check to see if each qinit file is present and then write the data
            for tfile in self.qinitfiles:
                try:
                    fname = "'%s'" % os.path.abspath(tfile[-1])
                except:
                    raise Warning("File %s was not found." % fname)
                    # raise MissingFile("file not found")
                self._out_file.write("\n%s  \n" % fname)
                self._out_file.write("%3i %3i \n" % tuple(tfile[:-1]))
        elif self.qinit_type >= 5 and self.qinit_type <= 9:
            self.data_write('epsilon')
            self.data_write("init_location")
            self.data_write("wave_family")
            self.data_write("angle")
            self.data_write("sigma")
        else:
            raise ValueError("Invalid qinit_type parameter %s." % self.qinit_type)
        self.close_data_file()