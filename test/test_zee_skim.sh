#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh

python3 src/main.py skim \
  --filepath data/DT_2024/EGamma/Run2024H.txt \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024H_Golden.json \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --out out_skim \
  --channel zee \
  --nThreads 8 \
  --nsteps 4 \
  --step 0 \
  --progress_bar \
  --correction_json data/corrections/summer24_corrections.json \
  --correction_key Run2024H

RET=$?

if [ $RET -ne 0 ]; then
  exit $RET
fi

python3 src/main.py hist \
  -hconf data/histograms/JECs.ini,data/histograms/kinematics.ini \
  -fl out_skim/J4PSkim_zee_0.root \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --channel zee \
  -loc \
  -pbar \
  --nThreads 8 \
  --out out_hist \
  --run_tag zee

RET=$?

exit $RET
