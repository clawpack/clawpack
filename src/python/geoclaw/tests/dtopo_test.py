#!/usr/bin/env python

r"""Test suite for dtopotools"""

import os
import urllib
import tarfile
import unittest
import tempfile
import glob

import clawpack.geoclaw.dtopotools as dtopotools

def get_bathy(url, destination=os.getcwd(), force=False):
    r"""Get bathymetry file located at `url`

    Will check downloaded file's suffix to see if the file needs to be extracted
    """

    file_name = os.path.basename(url)
    output_path = os.path.join(destination, file_name)
    if not os.path.exists(output_path) or force:
        print "Downloading %s to %s..." % (url, output_path)
        urllib.urlretrieve(url, output_path)
        print "Finished downloading."
    else:
        print "Skipping %s, file already exists." % file_name

    if tarfile.is_tarfile(output_path):
        with tarfile.open(output_path, mode="r:*") as tar_file:
            tar_file.extractall(path=destination)



class UCSBTest(unittest.TestCase):

    r""""""

    dtopo_file_check_path = ""

    def setUp(self):
        # Make temp directory for storage of input and output
        self.temp_path = tempfile.mkdtemp()

        # Fetch the UCSB fault specification if needed
        base_url = "http://users.ices.utexas.edu/~kyle/bathy/"
        get_bathy(os.path.join(base_url, 'usgs_honshu.subfault.tar.bz2'),
                  destination=self.temp_path)


    def test_download(self):
        
        assert os.path.exists(os.path.join(self.temp_path,"usgs_honshu.subfault")), \
               "Unable to download subfault specification file."

    def test_reading(self):
        pass


    def test_writing(self):
        pass


    def test_plot_subfaults(self):
        pass


    def test_plot_seafloor_deformation(self):
        pass


    def tearDown(self):
        paths = glob.glob(os.path.join(self.temp_path,"*"))
        for path in paths:
            os.remove(path)
        os.rmdir(self.temp_path)
        


class TestSingleSubfault(unittest.TestCase):

    r""""""

    def setUp(self):
        # Make temp directory for storage of input and output
        self.temp_path = tempfile.mkdtemp()

        self.test_subfault = dtopotools.SubFault()
        self.test_subfault.coordinate_specification = "centroid"
        self.test_subfault.dimensions = [1.0, 0.8]
        self.test_subfault.units = {"coordinates":"km", "dimensions":"km",
                                    "slip":"m", "depth":"km", "mu":"dyne/cm^2"}
        self.test_subfault.t = numpy.linspace(0.0, 100.0, 5)

        # subfault.faulttype = 'static'
        # subfault.dtopotype = 3
        # subfault.mx = 121
        # subfault.my = 81
        # subfault.xlower = -0.4
        # subfault.xupper = 0.6
        # subfault.ylower = -0.4
        # subfault.yupper = 0.4
        # subfault.t0 = 0.
        # subfault.tfinal = 100.
        # subfault.ntimes = 5

    def test_reading(self):
        subfault = dtopotools.SubFault(path=os.path.join(self.temp_path, 'test.subfault'))

        assert subfault == test_subfault


    def test_creation(self):
        
        subfault = dtopotools.Subfault(path)


    def test_compute_seafloor(self):
        pass

    
    def test_write(self):
        pass


    def test_creation(self):

        assert True, "Not implemented yet."

    def tearDown(self):
        # Remove contents of scratch directory
        paths = glob.glob(os.path.join(self.temp_path,"*"))
        for path in paths:
            os.remove(path)
        os.rmdir(self.temp_path)


class TestFaultPlotting(unittest.TestCase):

    r""""""

    # Skip this test, need a display
    __test__ = False

    def test_plot(self):
        assert False, "uhoh"


if __name__ == "__main__":
    unittest.main()

# def test_dtopo_plotting():
#     r"""Show plots of dtopo classes"""
#     assert False, "Uh oh..."

# def test_tohoku_ucsb():
#     """
#     Test reading a subfault file, plotting the slip distribution,
#     creating a dtopo file, reading the dtopo file back in, and
#     plotting the seafloor deformation.

#     The subfaults model the 11 March 2011 Great Tohoku earthquake.
#     """

#     figure(figsize=(15,9))

#     subplot(121)
#     plot_subfaults_tohoku_ucsb()    # read and plot subfaults
#     title("Slip on subfaults")

#     make_dtopo_tohoku_ucsb()        # make dtopo file of seafloor deformation

#     subplot(122)
#     plot_tohoku_ucsb()              # read dtopo file and plot deformation

#     fname = 'tohoku_ucsb.png'
#     savefig(fname)
#     print "Created ",fname

# def plot_subfaults_tohoku_ucsb():
#     """
#     Read subfaults for the UCSB model of the 11 March 2011 Tohoku event
#     and plot the slip distribution.
#     This subfault model came from the SUBFAULT FORMAT link at the bottom
#     of the page
#         http://www.geol.ucsb.edu/faculty/ji/big_earthquakes/2011/03/0311_v3/Honshu.html
#     and was developed by the UCSB group of Chen Ji, et al.
#     """

#     subfaults = dtopotools.read_subfault_model_ucsb('tohoku-ucsb.txt')
#     dtopotools.plot_subfaults(subfaults,slip_color=True,plot_rake=True)
#     return subfaults
    

    
# def make_dtopo_tohoku_ucsb():
#     """
#     Create dtopo data file for deformation of sea floor due to 
#     11 March 2011 Tohoku event, using the UCSB model.
#     """

#     subfault_fname = 'tohoku-ucsb.txt'
#     subfaults = dtopotools.read_subfault_model_ucsb(subfault_fname)

#     dtopo_fname = 'tohoku-ucsb.tt3'

#     if os.path.exists(dtopo_fname):
#         print "*** Not regenerating dtopo file (already exists): %s" \
#                 % dtopo_fname
#     else:
#         print "Using Okada model to create %s " % dtopo_fname
        
#         # Needed for extent of dtopo file:
#         xlower = 140.
#         xupper = 146.
#         ylower = 35.
#         yupper = 41.
        
#         # dtopo parameters for 1 min resolution:
#         mx = int((xupper - xlower)*60 + 1)
#         my = int((yupper - ylower)*60 + 1)
        
#         # Create dtopo_params dictionary with parameters for dtopo file: 
#         dtopo_params = {}
#         dtopo_params['fname'] = dtopo_fname
#         dtopo_params['faulttype'] = 'static'
#         dtopo_params['dtopotype'] = 3
#         dtopo_params['mx'] = mx
#         dtopo_params['my'] = my
#         dtopo_params['xlower'] = xlower
#         dtopo_params['xupper'] = xupper
#         dtopo_params['ylower'] = ylower
#         dtopo_params['yupper'] = yupper
#         dtopo_params['t0'] = 0.
#         dtopo_params['tfinal'] = 1.
#         dtopo_params['ntimes'] = 2
        
#         dtopo = dtopotools.make_dtopo_from_subfaults(subfaults, dtopo_params)
#         return dtopo

# def plot_tohoku_ucsb():
#     dtopo_fname = 'tohoku-ucsb.tt3'
#     dtopo = dtopotools.read_dtopo(dtopo_fname, 3)
#     x = dtopo.x
#     y = dtopo.y
#     dz_final = dtopo.dz_list[-1]
#     dtopotools.plot_dz_colors(x,y,dz_final,cmax_dz=16,dz_interval=1)
#     return dtopo


# #==================================================================

# def test_dtopo1():
#     """
#     Test reading a subfault file, plotting the slip distribution,
#     creating a dtopo file, reading the dtopo file back in, and
#     plotting the seafloor deformation.

#     Simple model with a single subfault.  
#     Vary the parameters to explore the Okada model.
#     """

#     figure(2,figsize=(15,9))

#     subplot(121)
#     subfaults = read_subfaults_dtopo1()    # read and plot subfaults
#     title("Slip on subfaults")
#     xlim(-0.4,0.4)
#     ylim(-0.4,0.4)

#     dtopo = make_dtopo_dtopo1()   # make dtopo file of seafloor deformation

#     subplot(122)
#     plot_dtopo1()              # read dtopo file and plot deformation

#     fname = 'dtopo1.png'
#     savefig(fname)
#     print "Created ",fname


# def read_subfaults_dtopo1():
#     """
#     Test data
#     """
#     fname_subfaults = 'dtopo1.csv'

#     # Format of subfault file:
#     columns = """longitude latitude depth length width strike dip rake slip
#        """.split()
#     defaults = {'latlong_location': 'top center'}
#     units = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}

#     subfaults = dtopotools.read_subfault_model(fname_subfaults, \
#                         columns=columns, units=units, \
#                         defaults = defaults, skiprows=1, delimiter=',')
                    
#     dtopotools.plot_subfaults(subfaults,slip_color=True,plot_rake=True)
        
#     return subfaults

# def make_dtopo_dtopo1():
#     """
#     Test data.
#     """

#     subfaults = read_subfaults_dtopo1()

#     dtopo_fname = 'dtopo1.tt3'

#     if os.path.exists(dtopo_fname):
#         print "*** Not regenerating dtopo file (already exists): %s" \
#                 % dtopo_fname
#     else:
#         print "Using Okada model to create %s " % dtopo_fname

#         # Needed for extent of dtopo file:
#         xlower = -0.4
#         xupper = 0.6
#         ylower = -0.4
#         yupper = 0.4

#         # number of grid points in dtopo file
#         mx = 121
#         my = 81

#         # Create dtopo_params dictionary with parameters for dtopo file: 
#         dtopo_params = {}
#         dtopo_params['fname'] = dtopo_fname
#         dtopo_params['faulttype'] = 'static'
#         dtopo_params['dtopotype'] = 3
#         dtopo_params['mx'] = mx
#         dtopo_params['my'] = my
#         dtopo_params['xlower'] = xlower
#         dtopo_params['xupper'] = xupper
#         dtopo_params['ylower'] = ylower
#         dtopo_params['yupper'] = yupper
#         dtopo_params['t0'] = 0.
#         dtopo_params['tfinal'] = 100.
#         dtopo_params['ntimes'] = 5

#         dtopo = dtopotools.make_dtopo_from_subfaults(subfaults, dtopo_params)

#         return dtopo

# def plot_dtopo1():
#     dtopo_fname = 'dtopo1.tt3'
#     dtopo = dtopotools.read_dtopo(dtopo_fname, 3)
#     x = dtopo.x
#     y = dtopo.y
#     dz_final = dtopo.dz_list[-1]
#     dtopotools.plot_dz_colors(x,y,dz_final,cmax_dz=16,dz_interval=1)
#     return dtopo

# if __name__=='__main__':

#     test_tohoku_ucsb()
#     test_dtopo1()
