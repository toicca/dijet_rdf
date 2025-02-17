#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh

python3 src/main.py skim \
  --filepath data/DT_2024/Muon/Run2024H.txt \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024H_Golden.json \
  --triggerpath data/triggerlists/ZJET_triggers_skim.txt \
  --out out_skim \
  --channel zjet \
  --nThreads 8 \
  --nsteps 20 \
  --step 1 \
  --correction_json data/corrections/summer24_corrections.json \
  --correction_key Run2024H \
  --progress_bar
