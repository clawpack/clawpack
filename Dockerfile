
ARG SOURCE_IMAGE=$SOURCE_IMAGE
FROM $SOURCE_IMAGE

COPY . /clawpack

## install Clawpack
ENV CLAW=/clawpack
ENV NETCDF4_DIR=/opt/conda
ENV FC=gfortran
ENV MPLBACKEND=Agg
# to make sure cython flags get appended rather than overwritten
ENV NPY_DISTUTILS_APPEND_FLAGS=1

# this is needed to find libraries when building geoclaw (particularly lapack)
ENV LIB_PATHS=/opt/conda/lib

WORKDIR /clawpack

# install clawpack
RUN pip install -e .

# install pytides
RUN pip install git+https://github.com/maritimeplanning/pytides.git@master \
      --no-cache-dir

# install nose
RUN conda update -n base conda
RUN conda install -yc conda-forge nose

WORKDIR /

ENTRYPOINT ["tini", "--", "/usr/bin/prepare.sh"]
