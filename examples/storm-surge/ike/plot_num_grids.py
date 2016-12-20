#!/usr/bin/env python

from __future__ import absolute_import
import os
import sys
import glob
import datetime

import numpy

# Plot customization
import matplotlib
from six.moves import range

# Markers and line widths
matplotlib.rcParams['lines.linewidth'] = 2.0
matplotlib.rcParams['lines.markersize'] = 6
matplotlib.rcParams['lines.markersize'] = 8

# Font Sizes
matplotlib.rcParams['font.size'] = 16
matplotlib.rcParams['axes.labelsize'] = 16
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rcParams['ytick.labelsize'] = 16

# DPI of output images
matplotlib.rcParams['savefig.dpi'] = 300

import matplotlib.pyplot as plt

convert2rgbfloat = lambda rgb: [value / 256.0 for value in rgb]
days2seconds = lambda days: days * 60.0**2 * 24.0
seconds2days = lambda seconds: seconds / (60.0**2 * 24.0)

def set_day_ticks(new_ticks=[-3, -2, -1, 0]):
    plt.xticks(new_ticks, [str(tick) for tick in new_ticks])

def set_cell_ticks():
    # plt.ticklabel_format(style='sci')
    locs,labels = plt.yticks()
    labels = locs / 1e6
    plt.yticks(locs,labels)
    # plt.yticks(new_ticks, [str(tick) for tick in new_ticks])

if __name__ == "__main__":

    output_path = "./_output"    
    if len(sys.argv) > 1:
        output_path = sys.argv[1]

    num_levels = 7
    # ADCIRC_num_nodes = 3331560
    landfall = datetime.datetime(2008,9,13 - 1,7) - datetime.datetime(2008,1,1,0)
    landfall = days2seconds(landfall.days) + landfall.seconds

    file_list = glob.glob(os.path.join(output_path,"fort.q*"))

    time = numpy.empty(len(file_list), dtype=float)
    num_grids = numpy.zeros((time.shape[0], num_levels), dtype=int)
    num_cells = numpy.zeros((time.shape[0], num_levels), dtype=int)

    for (n,path) in enumerate(file_list):
        # Read t file
        t_path = path[:-5] + "t" + path[-4:]
        t_file = open(t_path, 'r')
        time[n] = seconds2days(float(t_file.readline().split()[0]) - landfall)
        t_file.readline()
        t_file_num_grids = int(t_file.readline().split()[0])
        t_file.close()

        # Read q_file
        q_file = open(path, 'r')
        line = "\n"
        while line != "":
            line = q_file.readline()
            if "grid_number" in line:
                # print "grid number:", int(line.split()[0])
                level = int(q_file.readline().split()[0])
                num_grids[n,level - 1] += 1 
                mx = int(q_file.readline().split()[0])
                my = int(q_file.readline().split()[0])
                num_cells[n,level - 1] += mx * my

        q_file.close()

        # File checking
        if numpy.sum(num_grids[n,:]) != t_file_num_grids:
            raise Exception("Number of grids in fort.t* file and fort.q* file do not match.")

    # Plot cascading time histories per level
    colors = [ (value / 256.0, value / 256.0, value / 256.0) 
                                    for value in [247, 217, 189, 150, 115, 82, 37] ]
    proxy_artists = [plt.Rectangle((0, 0), 1, 1, fc=colors[level], 
            label="Level %s" % (str(level+1))) for level in range(num_levels)]

    # Number of grids
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.set_yscale('log')
    axes.stackplot(time ,num_grids.transpose(), colors=colors)
    axes.set_xlabel('Days from landfall')
    plt.subplots_adjust(left=0.13, bottom=0.12, right=0.90, top=0.90)
    axes.set_xlim([-2,0.75])
    axes.set_ylim([0,1e4])
    set_day_ticks()
    axes.set_ylabel('Number of Grids')
    axes.set_title("Number of Grids per Level in Time")
    axes.legend(proxy_artists, ["Level %s" % (str(level+1)) for level in range(num_levels)], loc=2)
    fig.savefig("num_grids.png")

    # Number of cells
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.set_yscale('log')
    axes.stackplot(time, num_cells.transpose(), colors=colors)
    axes.set_xlim([-2,0.75])
    axes.set_ylim([0,1e7])
    set_day_ticks()
    plt.subplots_adjust(left=0.13, bottom=0.12, right=0.90, top=0.90)
    axes.set_xlabel('Days from landfall')
    axes.set_ylabel('Number of Cells')
    axes.set_title("Number of Cells per Level in Time")
    axes.legend(proxy_artists, ["Level %s" % (str(level+1)) for level in range(num_levels)], loc=2)
    fig.savefig("num_cells.png")


    plt.show()


