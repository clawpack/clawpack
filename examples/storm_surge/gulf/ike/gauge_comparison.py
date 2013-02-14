#!/usr/bin/env python

import os
import sys
import datetime

import scipy.io
import matplotlib.pyplot as plt

import clawpack.clawutil.clawdata as clawdata
import clawpack.visclaw.gaugetools as gaugetools

# Conversion shortcuts
days2seconds = lambda days: days * 60.0**2 * 24.0
date2seconds = lambda date: days2seconds(date.days) + date.seconds
seconds2days = lambda secs: secs / (24.0 * 60.0**2)
min2deg = lambda minutes: minutes / 60.0
ft2m = lambda x:0.3048 * x

def read_Ike_gauge_data(base_path, skiprows=5, verbose=True):
    r"""Read the gauge info data file.

    Returns a dictionary for each gauge in the table.
      Keys: 'location': (tuple), 'depth': float, 'gauge_no': int
            'mean_water': (ndarray), 't': (ndarray)

    """
    stations = {}
    station_info_file = open(os.path.join(base_path,'Ike_Gauges_web.txt'),'r')

    # Skip past header
    for i in xrange(skiprows):
        station_info_file.readline()

    # Read in each station
    for line in station_info_file:
        data_line = line.split()
        if data_line[6] == "OK":
            stations[data_line[0]] = {
                    'location':[float(data_line[4]) + min2deg(float(data_line[5])),
                                float(data_line[2]) + min2deg(float(data_line[3]))],
                         'depth':float(data_line[8]) + float(data_line[9]),
                      'gauge_no':0}
            if data_line[1] == '-':
                stations[data_line[0]]['gauge_no'] = ord(data_line[0])
            else:
                stations[data_line[0]]['gauge_no'] = int(data_line[1])
            if verbose:
                print "Station %s: %s" % (data_line[0],stations[data_line[0]])
            
            # Load and extract real station data
            data = scipy.io.loadmat(os.path.join(base_path,'result_%s.mat' % data_line[0]))
            stations[data_line[0]]['t'] = data['yd_processed'][0,:]
            stations[data_line[0]]['mean_water'] = data['mean_water'].transpose()[0,:]

    station_info_file.close()

    return stations


def load_geoclaw_gauge_data(only_gauges=None, base_path="_output", verbose=True):
    r"""Load all gauge data in gauge file at base_path/fort.gauge

    Returns a dictionary of GaugeSolution objects keyed by their gauge numbers.
    """

    gauges = {}

    # Read in gauge.data file
    gauge_info_file = clawdata.GaugeData()
    gauge_info_file.read(data_path=base_path,file_name='gauges.data')

    if only_gauges is None:
        gauge_list = gauge_info_file.gauge_numbers
    else:
        gauge_list = only_gauges

    # Read in each gauge solution
    for (i,gauge_no) in enumerate(gauge_list):
        gauge = gaugetools.GaugeSolution(gauge_no, location=gauge_info_file.gauges[i][1:3])
        gauge.read(output_path=base_path, file_name='fort.gauge')
        gauges[gauge_no] = gauge
        if verbose:
            print "Read in GeoClaw gauge %s" % gauge_no

    return gauges


def plot_comparison(gauge_path,geoclaw_path):

    # Parameters
    surface_offset = 0.27
    landfall = []
    landfall.append(datetime.datetime(2008,9,13 + 1,7) 
                                                - datetime.datetime(2008,1,1,0))
    landfall.append(datetime.datetime(2008,9,13 - 1,7) 
                                                - datetime.datetime(2008,1,1,0))

    # Load gauge data
    kennedy_gauges = read_Ike_gauge_data(gauge_path)
    gauge_list = [gauge['gauge_no'] for gauge in kennedy_gauges.itervalues()]
    geoclaw_gauges = load_geoclaw_gauge_data(only_gauges=gauge_list,
                                             base_path=geoclaw_path)

    # Plot each matching gauge
    fig = plt.figure(figsize=(16,10),dpi=80)
    fig.suptitle('Surface from Sea Level')
    index = 0
    for (name,kennedy_gauge) in kennedy_gauges.iteritems():
        geoclaw_gauge = geoclaw_gauges[kennedy_gauge['gauge_no']]
        index = index + 1
        axes = fig.add_subplot(2,len(kennedy_gauges)/2,index)

        axes.plot(kennedy_gauge['t'] - seconds2days(date2seconds(landfall[0])),
                  kennedy_gauge['mean_water'] + kennedy_gauge['depth'], 'k')
        axes.plot(seconds2days(geoclaw_gauge.t - date2seconds(landfall[1])),
                  geoclaw_gauge.q[3,:] + surface_offset, 'r')

        axes.set_xlabel('Landfall Day')
        axes.set_ylabel('Surface (m)')
        axes.set_title("Station %s" % name)
        axes.set_xticks([-2,-1,0,1,2])
        axes.set_xticklabels([-2,-1,0,1,2])
        axes.set_xlim([-2,2])
        axes.set_ylim([-1,5])
        axes.grid(True)

    return fig


if __name__ == '__main__':
    kennedy_gauge_path = './Ike_gauge_data'
    geoclaw_output_path = "./_output"
    if len(sys.argv) > 1:
        geoclaw_output_path = sys.argv[1]
        if len(sys.argv) > 2:
            kennedy_gauge_path = sys.argv[2]

    # Plot Andrew Kennedy's gauge data versus corresponding GeoClaw data
    figure = plot_comparison(kennedy_gauge_path,geoclaw_output_path)
    plt.show()