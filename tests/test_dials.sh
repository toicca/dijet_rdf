
#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh

# Produce a DIALS json for the given files
python3 -m src.jec4prompt.main produce_dials \
    -ws jetmet \
    --json_file test_dials.json \
    --filepaths data/DT_2024/JetMET/Run2024I_fib.txt \
    --nThreads 8 

# Run over the different channels with the produced DIALS json, DCS, and golden json
GOLDEN_JSON=/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024I_Golden.json
DCS_JSON=/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/DCSOnly_JSONS/dailyDCSOnlyJSON/Collisions24_13p6TeV_378981_386951_DCSOnly_TkPx.json
DIALS_JSON=test_dials.json

python3 -m src.jec4prompt.main --tag golden skim \
  --filepath data/DT_2024/EGamma/Run2024I_fib.txt \
  --golden_json $GOLDEN_JSON \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --out out_skim_jsons \
  --channel photonjet \
  --nThreads 8 \
  --nsteps 1 \
  --step 0 \
  --progress_bar 

python3 -m src.jec4prompt.main hist \
    -hconf data/histograms/JECs.ini,data/histograms/kinematics.ini \
    -regions data/histograms/regions.ini \
    -fl out_skim_jsons/J4PSkimgolden_photonjet_0 \
    --triggerfile data/triggerlists/triggers_summer24.json \
    --channel photonjet \
    -loc \
    -pbar \
    --nThreads 8 \
    --out out_hist_jsons \
    --run_tag golden

python3 -m src.jec4prompt.main --tag dcs skim \
  --filepath data/DT_2024/EGamma/Run2024I_fib.txt \
  --golden_json $DCS_JSON \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --out out_skim_jsons \
  --channel photonjet \
  --nThreads 8 \
  --nsteps 1 \
  --step 0 \
  --progress_bar

python3 -m src.jec4prompt.main hist \
    -hconf data/histograms/JECs.ini,data/histograms/kinematics.ini \
    -regions data/histograms/regions.ini \
    -fl out_skim_jsons/J4PSkimdcs_photonjet_0 \
    --triggerfile data/triggerlists/triggers_summer24.json \
    --channel photonjet \
    -loc \
    -pbar \
    --nThreads 8 \
    --out out_hist_jsons \
    --run_tag dcs

python3 -m src.jec4prompt.main --tag dials skim \
  --filepath data/DT_2024/EGamma/Run2024I_fib.txt \
  --golden_json $DIALS_JSON \
  --triggerfile data/triggerlists/triggers_summer24.json \
  --out out_skim_jsons \
  --channel photonjet \
  --nThreads 8 \
  --nsteps 1 \
  --step 0 \
  --progress_bar

python3 -m src.jec4prompt.main hist \
    -hconf data/histograms/JECs.ini,data/histograms/kinematics.ini \
    -regions data/histograms/regions.ini \
    -fl out_skim_jsons/J4PSkimdials_photonjet_0 \
    --triggerfile data/triggerlists/triggers_summer24.json \
    --channel photonjet \
    -loc \
    -pbar \
    --nThreads 8 \
    --out out_hist_jsons \
    --run_tag dials

# python3 -m src.jec4prompt.main hist \
# -hconf data/histograms/JECs.ini,data/histograms/kinematics.ini \
# -fl out_skim/J4PSkim_dijet_0.root \
# --triggerfile data/triggerlists/triggers_summer24.json \
# --channel dijet \
# -loc \
# -pbar \
# --nThreads 8 \
# --out out_hist \
# --run_tag dijet

# RET=$?

# exit $RET