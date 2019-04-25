
ARG SOURCE_IMAGE=$SOURCE_IMAGE
FROM $SOURCE_IMAGE

RUN apt-get update -qq && \
    apt-get install -y \
        liblapack-pic \
        liblapack-dev \
        libproj-dev \
        proj-data \
        proj-bin \
        libgeos-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY . /clawpack

## install Clawpack
ENV CLAW=/clawpack
ENV NETCDF4_DIR=/usr/local
ENV FC=gfortran
ENV MPLBACKEND=Agg

# need to change shell in order for source command to work
SHELL ["/bin/bash", "-c"]

WORKDIR /clawpack

## currently pinning rhg_compute_tools until latest worker image version has
## more current version
RUN source activate worker && \
  pip install --upgrade pip && \
  pip install -e . && \
  pip install --upgrade matplotlib yolk3k \
   pytides rhg_compute_tools>=0.1.6 --no-cache-dir

WORKDIR /

ENTRYPOINT ["/usr/local/bin/dumb-init", "/usr/bin/prepare.sh"]
