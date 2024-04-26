#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh

python3 src/main.py --filepath data/testfiles.txt --selection_only --triggerpath data/testtrigger.txt \
    --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json \
    --output_path out --run_id rerun41b --nThreads 8 --number_of_files 4 --progress_bar 