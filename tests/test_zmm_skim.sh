#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh

python3 src/main.py skim \
  --filepath data/DT_2024/Muon/Run2024H.txt \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024H_Golden.json \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --out out_skim \
  --channel zmm \
  --nThreads 8 \
  --nsteps 1 \
  --step 0 \
  --correction_json data/corrections/summer24_corrections.json \
  --correction_key Run2024H \
  --progress_bar 

RET=$?

if [ $RET -ne 0 ]; then
  exit $RET
f

python3 src/main.py hist \
  -hconf data/histograms/JECs.ini,data/histograms/kinematics.ini,data/histograms/validation.ini \
  -fl out_skim/J4PSkim_zmm_0.root \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --channel zmm \
  -loc \
  -pbar \
  --nThreads 8 \
  --out out_hist \
  --run_tag zmm

RET=$?

exit $RET
