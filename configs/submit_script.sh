#!/bin/bash

CONFIG_FILE=$1

source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh

python3 src/main.py --config $CONFIG_FILE

eos cp *.root /eos/user/n/ntoikka/dijet_rdf/out
