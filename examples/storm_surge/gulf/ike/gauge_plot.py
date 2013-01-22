#!/usr/bin/env python

import os

import scipy.io
import matplotlib.pyplot as plt

import clawpack.visclaw.gaugetools as gaugetools

min2deg = lambda minutes: minutes / 60.0
ft2m = lambda x:0.3048 * x


# Parameters
base_path = './Ike_gauge_data'
Ike_landfall = 256

# Load in station info
stations = {}
station_info_file = open(os.path.join(base_path,'Ike_Gauges_web.txt'),'r')

# Skip past header
station_info_file.readline()
station_info_file.readline()
station_info_file.readline()
station_info_file.readline()
station_info_file.readline()

# Read in data lines
for line in station_info_file:
    data_line = line.split()
    if data_line[5] == "OK":
        stations[data_line[0]] = {
                'location':[float(data_line[3]) + min2deg(float(data_line[4])),
                            float(data_line[1]) + min2deg(float(data_line[2]))],
                     'depth':float(data_line[7]) + float(data_line[8])}
                   # 'depth':ft2m(float(data_line[6]))}

        print "Station %s: %s" % (data_line[0],stations[data_line[0]])
station_info_file.close()

# Plot data
fig_Hs = plt.figure()
fig_Hs.suptitle("Surface from Sea Level")

index = 0
for (name,station_data) in stations.iteritems():
    index = index + 1

    # Load and extract real station data
    data = scipy.io.loadmat(os.path.join(base_path,'result_%s.mat' % name))
    days = data['yd_processed'][0,:] - Ike_landfall
    mean_water = data['mean_water'].transpose()[0,:]

    # Load GeoClaw gauge data
    gauge = gaugetools.GaugeSolution(index,location=station_data['location'])
    gauge.read(output_path='./_output',file_name='fort.gauge')

    # Plot 
    axes = fig_Hs.add_subplot(2,len(stations)/2,index)
    axes.plot(days,mean_water+station_data['depth'],'k')
    axes.plot(gauge.t - Ike_landfall,gauge.q[3,:],'.b')
    axes.set_xlabel('Landfall Day')
    axes.set_ylabel('Height (m)')
    axes.set_title("Station %s" % name)
    axes.set_xticks([-1,0,1,2,3])
    axes.set_xticklabels([-1,0,1,2,3])
    axes.set_xlim([-1,3])
    axes.set_ylim([-1.0,5.0])
    axes.grid(True)

plt.show()