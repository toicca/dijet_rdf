# dijet_rdf
JEC4PROMPT backend for processing NanoAOD input.

## Setup
For running on lxplus use the LCG image, ie.,
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_107a/x86_64-el9-gcc14-opt/setup.sh
```

## Running the analysis
There are currently several scripts included in the repository for tasks such as
- Skimming
- Producing histograms from skimmed files
- Finding a JSON file to corresponding input files
- Finding the first and last run for given input files

For how to use these see README.md in the src/ directory.