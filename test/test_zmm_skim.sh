#!/bin/bash
python3 src/main.py skim \
  --filepath data/DT_2024/Muon/Run2024H.txt \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024H_Golden.json \
  --triggerpath data/triggerlists/ZJET_triggers_skim.txt \
  --out out_skim \
  --dataset zjet \
  --nThreads 8 \
  --nsteps 20 \
  --step 1 \
  --progress_bar
