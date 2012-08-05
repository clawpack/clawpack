#!/usr/bin/env python

import os
import compare_gauges

# Base directory
base_path = os.path.expandvars("$DATA_PATH/ml_2d/plane_wave/")
output_path = os.path.join(base_path,"gauge_comparisons")
if not os.path.exists(output_path):
    os.mkdir(output_path)

# Comparisons
comparisons = [["ml_2d_ia0_ba0_output","ml_2d_ia45_ba45_output"],
               ["ml_2d_ia0_ba45_output","ml_2d_ia45_ba0_output"],
               ["ml_2d_ia22_ba0_output","ml_2d_ia45_ba22_output"]]
comparison_paths = []
for path in comparisons:
    comparison_paths.append([os.path.join(base_path,path[0]),os.path.join(base_path,path[1])])
              
for (i,data_paths) in enumerate(comparison_paths):
    figures = compare_gauges.plot_gauges('all',data_paths,plot_vars=[[6,7],[6,7]],
                               titles=["Top Surface","Bottom Surface"],
                               kwargs={'figsize':(10,4)})
    for (num,figure) in figures.iteritems():
        name = os.path.split(data_paths[0])[-1][:-7]
        figure.savefig(os.path.join(output_path,"%s_%s_gauge.pdf" % (name,num)))
