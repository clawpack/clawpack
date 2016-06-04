r"""
Execute nosetests in all subdirectories, to run a series of quick
regression tests.

Sends output and result/errors to separate files to simplify checking
results and looking for errors.
"""

import os
import glob

import numpy

import clawpack.clawutil.test
import clawpack.pyclaw.util

# Clean library files whenever this module is used
if os.environ.has_key("CLAW"):
    CLAW = os.environ["CLAW"]
else:
    raise ValueError("Need to set CLAW environment variable.")

for lib_path in [os.path.join(CLAW,"amrclaw","src","2d"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow","multilayer"),
                 os.path.join(CLAW,"geoclaw","src","2d","shallow","surge")]:
    for path in glob.glob(os.path.join(lib_path,"*.o")):
        os.remove(path)
    for path in glob.glob(os.path.join(lib_path,"*.mod")):
        os.remove(path)


class GeoClawRegressionTest(clawpack.clawutil.test.ClawpackRegressionTest):

    r"""Base GeoClaw regression test setup derived from ClawpackRegressionTest

    """

    __doc__ += clawpack.pyclaw.util.add_parent_doc(
                                  clawpack.clawutil.test.ClawpackRegressionTest)


    def build_executable(self, executable_name="xgeoclaw"):
        r"""Build executable by running `make .exe` in test directory.

        Moves the resulting executable to the temporary directory.


        """

        super(GeoClawRegressionTest, self).build_executable(
                                                executable_name=executable_name)


    def check_fgmax(self, save=False):
        r"""Basic test to assert fgmax equality
        Currently just records sum of fg.h and of fg.s.

        :Input:
         - *save* (bool) - If *True* will save the output from this test to 
           the file *regresion_data.txt*.  Default is *False*.
        """

        from clawpack.geoclaw import fgmax_tools

        fg = fgmax_tools.FGmaxGrid()
        fname = os.path.join(self.temp_path, 'fgmax1.txt')
        fg.read_input_data(fname)
        fg.read_output(outdir=self.temp_path)

        data_sum = numpy.array([fg.h.sum(), fg.s.sum()])

        # Get (and save) regression comparison data
        regression_data_file = os.path.join(self.test_path, "regression_data",
                "regression_data_fgmax.txt")
        if save:
            numpy.savetxt(regression_data_file, data_sum)
        regression_sum = numpy.loadtxt(regression_data_file)

        # Compare data
        tolerance = 1e-14
        assert numpy.allclose(data_sum, regression_sum, tolerance), \
                "\n data: %s, \n expected: %s" % (data_sum, regression_sum)