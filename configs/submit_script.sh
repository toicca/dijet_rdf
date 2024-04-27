#!/bin/bash

export EOS_MGM_URL=root://eosuser.cern.ch
CONFIG_FILE=$1

source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh
echo "Running on $(hostname)"

eos cp -r /eos/user/n/ntoikka/dijet_rdf .
cd dijet_rdf
python3 src/main.py --config $CONFIG_FILE > out.txt

cp *.root /eos/user/n/ntoikka/dijet_rdf/out/newest_out.root
cp out.txt /eos/user/n/ntoikka/dijet_rdf/out_2.txt
