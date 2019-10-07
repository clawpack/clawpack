
ARG SOURCE_IMAGE=$SOURCE_IMAGE
FROM $SOURCE_IMAGE

COPY . /clawpack

## install Clawpack
ENV CLAW=/clawpack
ENV NETCDF4_DIR=/opt/conda
ENV FC=gfortran
ENV MPLBACKEND=Agg

WORKDIR /clawpack

RUN pip install -e . && \
    pip install yolk3k nose \
      git+https://github.com/maritimeplanning/pytides.git@master \
      --no-cache-dir

WORKDIR /

ENTRYPOINT ["tini", "--", "/usr/bin/prepare.sh"]
