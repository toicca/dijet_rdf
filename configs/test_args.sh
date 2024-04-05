#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh

python3 src/main.py --filepath data/testfiles.txt --triggerpath data/testtrigger.txt --output_path out --run_id test_args --nThreads 8 --number_of_files 5 --progress_bar --is_local
