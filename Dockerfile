
FROM rhodium/worker:v0.2.5

COPY . /clawpack

## install Clawpack
ENV CLAW=/clawpack

# need to change shell in order for source command to work
SHELL ["/bin/bash", "-c"]

RUN source activate worker && \
  pip install -e .

ENTRYPOINT ["/usr/local/bin/dumb-init", "/usr/bin/prepare.sh"]
