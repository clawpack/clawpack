
ARG SOURCE_IMAGE=$SOURCE_IMAGE
FROM $SOURCE_IMAGE

RUN apt-get update -qq && \
    apt-get install -y \
        liblapack-pic \
        liblapack-dev \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY . /clawpack

## install Clawpack
ENV CLAW=/clawpack
ENV NETCDF4_DIR=/opt/conda
ENV FC=gfortran
ENV MPLBACKEND=Agg

WORKDIR /clawpack

RUN conda install -c conda-forge nose

RUN pip install -e . && \
    pip install yolk3k \
      git+https://github.com/maritimeplanning/pytides.git@master \
      --no-cache-dir

WORKDIR /

ENTRYPOINT ["tini", "--", "/usr/bin/prepare.sh"]
