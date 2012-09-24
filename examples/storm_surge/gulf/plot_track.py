#!/usr/bin/env python

r"""Plot hurricane track given a storm.data file"""

import numpy as np

import matplotlib.pyplot as plt

sec2days = lambda seconds:seconds / (60.0**2 * 24.0)

colors = ['r','b']

def plot_track(t,x,y,wind_radius,wind_speed,Pc):

    divide = (np.max(Pc) + np.min(Pc)) / 2.0

    fig = plt.figure(1)
    axes = fig.add_subplot(111)
    indices = Pc < divide
    axes.scatter(x[indices],y[indices],color='r',marker='o')
    indices = Pc >= divide
    axes.scatter(x[indices],y[indices],color='b',marker='o')
    axes.set_title("Track - Hurricane Ike")

    fig = plt.figure(2,figsize=(8*3,6))
    axes = fig.add_subplot(131)
    axes.plot(sec2days(t),wind_speed)
    axes.set_title("Maximum Wind Speed - Hurricane Ike")

    axes = fig.add_subplot(132)
    axes.plot(sec2days(t),wind_radius)
    axes.set_title("Maximum Wind Radius - Hurricane Ike")

    axes = fig.add_subplot(133)
    axes.plot(sec2days(t),Pc)
    axes.plot(sec2days(t),np.ones(t.shape) * divide,'k--')
    axes.set_title("Central Pressure - Hurricane Ike")



if __name__ == "__main__":
    data = np.loadtxt('_output/storm.data')
    plot_track(data[:,0],data[:,1],data[:,2],data[:,-3],data[:,-2],data[:,-1])

    plt.show()