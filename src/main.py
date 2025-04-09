from cmsdials.auth.bearer import Credentials
from cmsdials import Dials
from cmsdials.filters import (
    RunFilters,
    MLModelsIndexFilters,
    MLBadLumisectionFilters
)

# If you know the exact dataset you are looking
# it is even better, you can use RunFilters(dataset="...")
# otherwise, use the regex.
DATASET__REGEX = "JetMET./Run2024.*-PromptReco.*/DQMIO"
RUNS = [383811, 383812, 383813, 383814, 383830]

# Authenticate
creds = Credentials.from_creds_file()
dials = Dials(creds, workspace="jetmet")

# Since your starting point is a NanoAOD file
# DAS: dataset status=* dataset=/JetMET0/Run2024*-PromptReco*/NANOAOD
# for each dataset -> many files -> for each file -> many runs
runs_in_dials = dials.run.list_all(RunFilters(page_size=500, dataset__regex=DATASET__REGEX))
# runs_in_dials = runs_in_dials.to_pandas()
# print(runs_in_dials.results)
runs =  [run.run_number for run in runs_in_dials.results]
print(runs)
print(len(runs))
exit()

# You possibly have list of runs for the NANOAOD file you are analyzing
# Beaware that DIALS may not have ingested these runs yet, so check the filtered object
# to be sure if something is missing.
runs_to_analyze = runs_in_dials[runs_in_dials.run_number.isin(RUNS)]

# At this point, you have: dataset ids and runs available in dials
# So you just need to know the models you want to use
#
# Filtering active=True means: Models that are currently active and predicting new data
#
# Right now we have limited metadata for filtering, we are working in enhancing the ML part!
models_in_dials = dials.ml_models_index.list_all(MLModelsIndexFilters(page_size=500, active=True))
models_in_dials = models_in_dials.to_pandas()

# MODE 1: Requesting raw predictions
raw_preds = dials.ml_bad_lumis.list_all(
    MLBadLumisectionFilters(
        page_size=500,
        model_id__in=models_in_dials.model_id.tolist(),
        dataset_id__in=runs_to_analyze.dataset_id.tolist(),
        run_number__in=runs_to_analyze.run_number.tolist()
    )
)
raw_preds = raw_preds.to_pandas()
print(raw_preds)

# MODE 2: Requesting predictions in the certification format
cert_json = dials.ml_bad_lumis.cert_json(
    model_id__in=models_in_dials.model_id.tolist(),
    dataset_id__in=runs_to_analyze.dataset_id.tolist(),
    run_number__in=runs_to_analyze.run_number.tolist()
)
print(cert_json)

# MODE 3: Requesting predictions in the golden-json-like format
mlgolden_json = dials.ml_bad_lumis.golden_json(
    model_id__in=models_in_dials.model_id.tolist(),
    dataset_id__in=runs_to_analyze.dataset_id.tolist(),
    run_number__in=runs_to_analyze.run_number.tolist()
)
print(mlgolden_json)

import json
with open("mlgolden.json", "w") as f:
    json.dump(mlgolden_json, f)
