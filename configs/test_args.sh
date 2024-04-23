#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_105a/x86_64-el9-gcc13-opt/setup.sh

python3 src/main.py --filepath data/testfiles.txt --triggerpath data/testtrigger.txt --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/DCSOnly_JSONS/dailyDCSOnlyJSON/Collisions24_13p6TeV_378981_379530_DCSOnly_TkPx.json --output_path out --run_id rerun32 --nThreads 8 --number_of_files 2 --progress_bar 