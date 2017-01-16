#!/bin/bash
env
export LDFLAGS="$LDFLAGS -Wl,-headerpad_max_install_names -undefined dynamic_lookup -bundle"
$PYTHON setup.py install

