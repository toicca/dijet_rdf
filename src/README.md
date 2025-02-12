# Running JEC4PROMPT
All scripts are ran through `main.py`, which has command line parser for each of the scripts provided.

## Skimming
Skims are created with the command `python3 main.py skim`.

The command requires
- input files passed with the flag `--filelist` or `--filepaths`. `--filelists` takes as input a commma separated list of root files, ie. `--filelist file1.root,file2.root`, and `--filepaths` takes as input a `.txt` file with a list of root files, ie. `--filepaths data/DT_2024/JetMET/Run2024H.txt`. If the files are located on the same filesystem as where the commands are run, the flag `--is_local` needs to be passed to the script.
- output directory passed with the `--out`, ie. `--out out_skim/`
- channel that the skim should be based on, ie. `--dataset dijet`.

Additionally, you can specify
- `--progress_bar` to follow the progress and performance of the skimming.
- `--triggerpath` or `--triggerlist` to skim based on triggers. The flags are labeled similar to `--filelist` and `--filepaths`, ie. you can use either `--triggerpath data/triggerlists/EGM_triggers_skim.txt` or `--triggerlist HLT_ZeroBias,HLT_PFJet40`
- `--nThreads` to specify number of threads for processing, ie. `--nThreads 8`
- `--golden_json` to specify a JSON to filter the runs and lumisections, ie. `--golden_json /eos/user/c/cmsdqm/www/CAF/certification/Collisions24/2024I_Golden.json`
- `--run_range` to specify the runs when the data was collected and to calculate the luminosity collected in that range, ie. `--run_range 384069,384128`. This option requires a JSON passed with the `--golden_json`.
- `--nsteps` and `--step` to split the input filelist to `--nsteps` and to process the step number `--step`, ie. `--nsteps 10 --step 0` will split the input files to 10 sets and process the first set.
- `--correction_json` and `--correction_key` to define the `data/` correction json to use and which corrections to use from the file, ie. `--correction_json data/corrections/summer24_corrections.json --correction_key 2024H`
