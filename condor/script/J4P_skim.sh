#!/bin/bash
export INFILE=$1
export OUTPATH=$2
export STEP=$3
export NSTEPS=$4
export CHANNEL=$5
export TRG=$6
export X509_USER_PROXY=$7

export EOS_MGM_URL="root://eosuser.cern.ch"

export EXTRA_CLING_ARGS='-O3'

source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh

eos cp -r $OUTPATH/data ./
eos cp -r $OUTPATH/src ./

python3 src/main.py skim \
  --filepath $INFILE \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024H_Golden.json \
  --triggerpath data/triggerlists/${TRG}_triggers_skim.txt \
  --out $OUTPATH/out_skim_condor \
  --channel $CHANNEL \
  --nThreads 4 \
  --nsteps $NSTEPS \
  --step $STEP \
  --correction_json data/corrections/summer24_corrections.json \
  --correction_key Run2024H

RET=$?

if [ $RET -ne 0 ]; then
  exit $RET
fi