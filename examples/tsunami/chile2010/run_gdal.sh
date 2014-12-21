#!/bin/tcsh

set file = frame000$1fig1
set ext = png

setenv GDAL_DATA /opt/local/share/gdal

gdalinfo $file.$ext

gdal_translate -of VRT \
    -a_srs EPSG:4326 \
    -gcp 0        0 -120   0      \
    -gcp 2400     0  -60   0      \
    -gcp 2400  2400  -60 -60    \
    -90 \
    $file.$ext "$file"_tmp.vrt

gdalwarp -of VRT \
       -t_srs EPSG:4326 \
       -overwrite \
	"$file"_tmp.vrt $file.vrt


gdal2tiles.py --zoom=2-6 --profile=geodetic  --force-kml --resampling=near $file.vrt
