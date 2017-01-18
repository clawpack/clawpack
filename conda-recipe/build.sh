#!/bin/bash
if [[ $(uname) == Darwin ]];
then
    export LDFLAGS="$LDFLAGS -Wl,-headerpad_max_install_names -undefined dynamic_lookup -bundle"
fi
$PYTHON setup.py install
