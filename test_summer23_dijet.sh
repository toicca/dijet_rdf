#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh

python3 src/main.py skim \
  --filepath data/DT_2023/JetMET/Run2023Cv123.txt \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json \
  --triggerfile data/triggerlists/triggers_summer23.json \
  --out out_skim_23Cv123 \
  --channel dijet \
  --nThreads 8 \
  --nsteps 1 \
  --step 0 \
  --correction_json data/corrections/summer23_corrections.json \
  --correction_key Run2023Cv123 \
  --progress_bar

python3 src/main.py skim \
  --filepath data/DT_2023/JetMET/Run2023Cv4.txt \
  --golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json \
  --triggerfile data/triggerlists/triggers_summer23.json \
  --out out_skim_23Cv4 \
  --channel dijet \
  --nThreads 8 \
  --nsteps 1 \
  --step 0 \
  --correction_json data/corrections/summer23_corrections.json \
  --correction_key Run2023Cv4 \
  --progress_bar

RET=$?

if [ $RET -ne 0 ]; then
  exit $RET
fi

python3 src/main.py hist \
  -hconf data/histograms/JECs.ini,data/histograms/kinematics.ini,data/histograms/validation.ini \
  -fl out_skim/J4PSkim_dijet_0.root \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --channel dijet \
  -loc \
  -pbar \
  --nThreads 8 \
  --out out_hist \
  --run_tag dijet

RET=$?

exit $RET