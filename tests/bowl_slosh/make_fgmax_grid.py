
"""
Create fgmax_grid.txt input file
"""

from __future__ import absolute_import
from clawpack.geoclaw import fgmax_tools
import os


def make_fgmax_grid1(datadir):
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2       # will specify a 2d grid of points
    fg.x1 = -2.
    fg.x2 = 2.
    fg.y1 = -2.
    fg.y2 = 2.
    fg.dx = 0.1
    fg.tstart_max = 0.        # when to start monitoring max values
    fg.tend_max = 1.e10       # when to stop monitoring max values
    fg.dt_check = 0.1         # target time (sec) increment between updating
                               # max values
    fg.min_level_check = 2    # which levels to monitor max on
    fg.arrival_tol = 1.e-2    # tolerance for flagging arrival

    fg.input_file_name = os.path.join(datadir, 'fgmax1.txt')
    fg.write_input_data()


if __name__ == "__main__":
    make_fgmax_grid1('.')


