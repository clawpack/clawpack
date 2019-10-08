
ARG SOURCE_IMAGE=$SOURCE_IMAGE
FROM $SOURCE_IMAGE

COPY . /clawpack

## install Clawpack
ENV CLAW=/clawpack
ENV NETCDF4_DIR=/opt/conda
ENV FC=gfortran
ENV MPLBACKEND=Agg
ENV LIB_PATHS=/opt/conda/lib

WORKDIR /clawpack

# need these to pass tests (and maybe we need them for our geoclaw runs but not
# sure)
RUN conda install -y -c conda-forge lapack nose h5py

WORKDIR /

ENTRYPOINT ["tini", "--", "/usr/bin/prepare.sh"]
