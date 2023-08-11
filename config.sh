#!/bin/bash

ulimit -Ss 64000
NEWLIMIT=$(ulimit -Ss)
echo "limit is set to $NEWLIMIT"

pip3 install -r requirements.txt
python3 setup.py build_ext --inplace
