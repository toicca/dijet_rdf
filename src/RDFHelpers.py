import ROOT
import json
import numpy as np
from typing import List
import argparse
import configparser

#====================================================
#
# Functions for files and parsing
#
#====================================================

def file_read_lines(file: str, find_ROOT: bool = False) -> List[str]:
    if find_ROOT:
        with open(file) as f:
            return [line.strip() for line in f.readlines() if line.strip().endswith(".root")]
    else:
        with open(file) as f:
            return [line.strip() for line in f.readlines()]
    
def read_config_file(config_file: str) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)
    return config
    
def parse_arguments():
    parser = argparse.ArgumentParser(description='JEC4PROMPT Analyzer')
    
    # General config
    filepath_group = parser.add_mutually_exclusive_group(required=True)
    filepath_group.add_argument('-conf', '--config', type=str, help='Path to the config file. If set, overrides all other options')
    filepath_group.add_argument('-fp', '--filepath', type=str, help='Path to the file list')
    filepath_group.add_argument('-fl', '--filelist', type=str, help='Input files separated by commas')
    trigger_group = parser.add_mutually_exclusive_group()
    trigger_group.add_argument('-tp', '--triggerpath', type=str, help='Path to the trigger list')
    trigger_group.add_argument('-tl','--triggerlist', type=str, help='Input files separated by commas')
    parser.add_argument('-nof', '--number_of_files', type=int, default=-1, help='How many files to be processed. -1 for all files in the list')
    parser.add_argument('-loc', '--is_local', action='store_true', help='Run locally. If not set will append root://cms-xrd-global.cern.ch/ to the start of file names')
    parser.add_argument('-out', '--output_path', type=str, help='Path where to write output files')
    parser.add_argument('-id', '--run_id', type=str, help='Run identifier such as date or version of the software included in output file names')
    parser.add_argument('-MC', '--is_MC', action='store_true', help='Set if running on MC')

    # Corrections and filtering files
    parser.add_argument('-gjson', '--golden_json', type=str, default='', help='Path to the golden JSON file') # good job son
    parser.add_argument('-jetvm', '--jetvetomap', type=str, help='Path to the jetvetomap file')
    parser.add_argument('-L1', '--L1FastJet', type=str, help='Path to the L1FastJet txt correction file (legacy option)')
    parser.add_argument('-L2Rel', '--L2Relative', type=str, help='Path to the L2Relative txt correction file')
    parser.add_argument('-L2Res', '--L2L3Residual', type=str, help='Path to the L2L3Residual txt correction file')
    parser.add_argument('-JER', '--JER', type=str, help='Path to the JER txt correction file')
    parser.add_argument('-JER_SF', '--JER_SF', type=str, help='Path to the JER scale factor txt correction file')

    # Performance and logging
    parser.add_argument('-nThreads', '--nThreads', type=int, default=2, help='Number of threads to use')
    parser.add_argument('-verb', '--verbosity', type=int, choices=[0,1,2], help='Verbosity level')   # TODO
    parser.add_argument('-pbar', '--progress_bar', action='store_true', help='Show progress bar')
    parser.add_argument('-cfrep', '--cutflow_report', action='store_true', help='Print cutflow report') # TODO

    # Parse command line arguments, overriding config file values
    args = parser.parse_args()
    
    # If config file is set, override all arguments with ones set in there
    if args.config:
        config = read_config_file(args.config)
        for arg, value in config["GENERAL"].items():
            # Do a type conversion for the option
            if arg == "number_of_files" or arg == "nThreads" or arg == "verbosity" or arg == "is_local" or arg == "is_mc" or arg == "progress_bar" or arg == "cutflow_report":
                if value == "":
                    setattr(args, arg, 0)
                else:
                    setattr(args, arg, int(value))
            else:
                setattr(args, arg, value)
                
    # Split the file list and trigger list if they are given as a string
    if args.filelist:
        args.filelist = args.filelist.split(",")
    elif args.filepath:
        args.filepath = file_read_lines(args.filepath, find_ROOT=True)
        
    if args.triggerlist:
        args.triggerlist = args.triggerlist.split(",")
    elif args.triggerpath:
        args.triggerpath = file_read_lines(args.triggerpath)

    return args


#====================================================
#
# Bins for histograms
#
#====================================================
def get_bins() -> dict:
    bins = {}
    
    # These pt and eta bins are some JEC bins
    bins["pt"] = {}
    bins["pt"]["bins"] = np.array((1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133,
                    153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737,
                    790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
                    2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
                    4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000), dtype=float)
    bins["pt"]["n"] = len(bins["pt"]["bins"]) - 1

    bins["eta"] = {}
    bins["eta"]["bins"] = np.array((-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314,
    -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566,
    -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522,
    -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609,
    0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74,
    1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839,
    4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191), dtype=float)
    bins["eta"]["n"] = len(bins["eta"]["bins"]) - 1
    
    bins["phi"] = {}
    bins["phi"]["bins"] = np.linspace(-np.pi, np.pi, 73, dtype=float)
    bins["phi"]["n"] = len(bins["phi"]["bins"]) - 1
    
    bins["mjj"] = {}
    bins["mjj"]["bins"] = np.linspace(200, 10000, 200, dtype=float)
    bins["mjj"]["n"] = len(bins["mjj"]["bins"]) - 1
    
    bins["deltaEta"] = {}
    bins["deltaEta"]["bins"] = np.linspace(0, 10, 100, dtype=float)
    bins["deltaEta"]["n"] = len(bins["deltaEta"]["bins"]) - 1
    
    bins["deltaR"] = {}
    bins["deltaR"]["bins"] = np.linspace(0, 10, 100, dtype=float)
    bins["deltaR"]["n"] = len(bins["deltaR"]["bins"]) - 1
    
    bins["response"] = {}
    bins["response"]["bins"] = np.linspace(0, 2, 100, dtype=float)
    bins["response"]["n"] = len(bins["response"]["bins"]) - 1
    
    bins["asymmetry"] = {}
    bins["asymmetry"]["bins"] = np.linspace(-1, 1, 100, dtype=float)
    bins["asymmetry"]["n"] = len(bins["asymmetry"]["bins"]) - 1
    
    bins["runs"] = {}
    bins["runs"]["bins"] = np.linspace(355065, 391370, int((391370-355065) / 10000), dtype=float)
    bins["runs"]["n"] = len(bins["runs"]["bins"]) - 1
    
    bins["bx"] = {}
    bins["bx"]["bins"] = np.linspace(0, 3564, 3564, dtype=float)
    bins["bx"]["n"] = len(bins["bx"]["bins"]) - 1
    
    bins["lumi"] = {}
    bins["lumi"]["bins"] = np.linspace(0, 1600, 1600, dtype=float)
    bins["lumi"]["n"] = len(bins["lumi"]["bins"]) - 1
    
    return bins

if __name__ == "__main__":
    # Test for the argument parser
    args = parse_arguments()
    print(args)
    if None:
        print("yeeters")