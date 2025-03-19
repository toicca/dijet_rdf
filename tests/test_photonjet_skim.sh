#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh

python3 -m src.jec4prompt.main skim \
  --filepath data/DT_2024/EGamma/Run2024H.txt \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024H_Golden.json \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --out out_skim \
  --channel photonjet \
  --nThreads 8 \
  --nsteps 10 \
  --step 0 \
  --correction_json data/corrections/summer24_corrections.json \
  --correction_key Run2024H \
  --progress_bar

RET=$?

if [ $RET -ne 0 ]; then
  exit $RET
fi

python3 -m src.jec4prompt.main hist \
  -hconf data/histograms/kinematics.ini,data/histograms/JECs.ini \
  -fl out_skim/J4PSkim_photonjet_0.root \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --channel photonjet \
  -loc \
  -pbar \
  --nThreads 8 \
  --out out_hist \
  --run_tag photonjet

RET=$?

exit $RET