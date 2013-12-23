
6 nested grids to debug new rectintegral.f90

Note that recursive Makefiles don't work with -f flag for .output, so need
to do the following:

make .exe -f Makefile_new
make data
make output -f Makefile_new
