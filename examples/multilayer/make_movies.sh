#!/bin/sh

# Make all movies
delay=10
loop=1

for simulation in ml_2d_ia0_ba0 ml_2d_ia0_ba45 ml_2d_ia22_ba0 ml_2d_ia45_ba0 ml_2d_ia45_ba22 ml_2d_ia45_ba45 ; do
    for figure in 0 1 ; do
        convert -delay ${delay} -loop ${loop} ${DATA_PATH}/ml_2d/plane_wave/${simulation}_plots/frame*fig${figure}.png ${DATA_PATH}/ml_2d/plane_wave/movies/${simulation}_fig${figure}.gif
    done
done
