# dijet_rdf
JEC4PROMPT backend for processing NanoAOD input.

## Setup
For running on lxplus use the LCG image, ie.,
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh
```

The toolkit can be installed (in editable mode) with
```
pip install -e .
```
or run directly with
```
python3 -m src.jec4prompt.main
```

For faster processing with RDF you can also include
```
export EXTRA_CLING_ARGS='-O3'
```

## Analysis metadata and inputs
Information about input files, corrections, trigger lists, etc., are hosted in the `data/` directory.

### Input files
Input files are listed in `.txt` files in directories such as `data/MC_Summer24` and `data/DT_2024`. The text files work directly as input for the skimming script.

### Corrections
It is recommended to use the JSONs provided by the JSONPOG hosted at `/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/` for corrections. Corrections are then defined through these files and the stacks in those files with `.jsons` such as `data/corrections/summer23_corrections.json`, and the correction files themselves are taken from `/cvmfs`.

### Trigger lists
When skimming a large part of the size reduction comes from requiring triggers to trigger. Triggers that should be used in a skim are defined in `.json` files in `data/triggerlists/`. The keys in the files define what triggers are used in _skims_, and the `cut` key defines further cuts per trigger when eg. creating histograms. The string given to `cut` in a trigger is directly used as an argument in an RDFs `.Filter()`, which means that it is possible to combine boolean expressions with the trigger to define eg. regions where the trigger is efficient / doesn't overlap with other analysis triggers.

## Running the analysis
There are currently several scripts included in the repository for tasks such as
- Skimming
- Producing histograms from skimmed files
- Finding a JSON file to corresponding input files
- Finding the first and last run for given input files

For how to use these see README.md in the src/ directory.