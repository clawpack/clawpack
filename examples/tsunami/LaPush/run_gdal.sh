#!/bin/tcsh

set file = frame0012fig1
set ext = png

setenv GDAL_DATA /opt/local/share/gdal

gdalinfo $file.$ext

gdal_translate -of VRT \
    -a_srs EPSG:4326 \
    -gcp 0 0      -124.68 47.96  \
    -gcp 9920 0   -124.55 47.96 \
    -gcp  9920 7630 -124.55 47.86 \
    -90 $file.$ext $file.vrt

 gdalwarp -of VRT \
     -t_srs EPSG:4326 \
     $file.vrt "$file".vrt

gdal2tiles.py --profile=geodetic  --force-kml --resampling=near "$file".vrt
