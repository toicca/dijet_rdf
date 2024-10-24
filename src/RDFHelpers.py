import ROOT
import json
import numpy as np
from typing import List, Dict
import argparse
import configparser

#====================================================
#
# Functions for files and parsing
#
#====================================================
def find_era(files: List[str]) -> str:
    eras = {}
    for file in files:
        idx = file.find("Run")
        if idx != -1:
            era = file[idx:idx+8]
            if era in eras:
                eras[era] += 1
            else:
                eras[era] = 1
    print(f"Found eras {eras} from files")
    print(f"Most common era: {max(eras, key=eras.get)}")
    return max(eras, key=eras.get)

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

def read_trigger_config(config_file: str) -> Dict:
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)

    triggers = {}
    for section in config.sections():
        if "filter" not in config[section]:
            triggers[section] = section
        else:
            triggers[section] = "(" + config[section]["filter"] + " && " + section + ")"

    return triggers
    
def read_correction_config(config_file: str) -> Dict:
    types_in_ROOT = {
        "RVec<float>": "ROOT::VecOps::RVec<float>",
    }

    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)

    corrections = {}
    for section in config.sections():
        corrections[section] = {}
        corrections[section]["file"] = config[section]["file"]

        types = config[section]["types"].split(",")
        types = [types_in_ROOT[type] if type in types_in_ROOT else type for type in types]
        variables = config[section]["variables"].split(",")

        func_call = f"ROOT::VecOps::RVec<float> getJEC{section}("
        for t, v in zip(types, variables):
            func_call += f"{t} {v}, "
        func_call = func_call[:-2] + ")"

        corrections[section]["func_call"] = func_call

        rdf_call = f"getJEC{section}("
        for v in variables:
            rdf_call += f"{v}, "
        rdf_call = rdf_call[:-2] + ")"
        corrections[section]["rdf_call"] = rdf_call

        eval_call = "{"
        for t, v in zip(types, variables):
            if "vec" in t.lower():
                eval_call += v + "[i], "
            else:
                eval_call += v + ", "
        eval_call = eval_call[:-2] + "}"
        corrections[section]["eval_call"] = eval_call

    return corrections

def parse_arguments():
    parser = argparse.ArgumentParser(description='JEC4PROMPT Toolkit: \
            https://github.com/toicca/dijet_rdf/tree/main')
    
    subparsers = parser.add_subparsers(dest='subparser_name')

    # Histogram config: Histogram producer from skimmed files
    hist_parser = subparsers.add_parser('hist', help='Produce histograms from skimmed files.')
    hist_files = hist_parser.add_mutually_exclusive_group(required=True)
    hist_files.add_argument('-c', '--config', type=str, help='Path to the config file. If set, \
            overrides all other options.')
    hist_files.add_argument('-fp', '--filepath', type=str, help='Path to the file list.')
    hist_files.add_argument('-fl', '--filelist', type=str, help='Input files separated by commas.')
    hist_parser.add_argument('-loc', '--is_local', action='store_true', help='Run locally. If not \
            set will append root://cms-xrd-global.cern.ch/ \
            to the start of file names.')
    hist_parser.add_argument('-pbar', '--progress_bar', action='store_true',
            help='Show progress bar.')
    hist_parser.add_argument('-hconf', '--hist_config', type=str, help='Path to the histogram \
            config file.')
    hist_parser.add_argument("--run_range", type=str, help="Run range of the given input files \
            (run_min and run_max separated by a comma)")
    hist_parser.add_argument("--run_tag", type=str, help="Run tag")
    hist_parser.add_argument("--nThreads", type=int, help="Number of threads to be used \
            for multithreading")
    hist_parser.add_argument("--out", type=str, required=True, default="", help="Output path \
            (output file name included)")
    
    # Find JSON config
    find_json_parser = subparsers.add_parser('find_json', help='Find JSON File appropriate for given run')
    find_json_parser.add_argument("--json_files", required=True, type=str, help="Comma separated \
            list of json files")
    find_json_parser.add_argument("--run", required=True, type=str, help="Run number")
    find_json_parser.add_argument("--out", required=False, type=str, help="Output file")

    # Find newest config
    find_newest_parser = subparsers.add_parser('find_newest', help="Find newest output file in \
            the subdirectories of given root directory")
    find_newest_parser.add_argument("--root_directory", type=str, help="Directory to search for \
            files in")
    find_newest_parser.add_argument("--starts_with", type=str, help="Choose a prefix for the files \
            to search for")
    find_newest_parser.add_argument("--ends_with", type=str, help="Choose a suffix for the files \
            to search for")
    find_newest_parser.add_argument("--spaces", action="store_true", help="Use spaces instead of \
            commas to separate the file paths")
    find_newest_parser.add_argument("--max_depth", type=int, help="Depth of files to search for \
            in the directory tree (default: None)")

    # Find range config
    find_range_parser = subparsers.add_parser('find_range', help="Find run range of \
                                                given input files")
    find_range_files = find_range_parser.add_mutually_exclusive_group(required=True)
    find_range_files.add_argument("--filelist", type=str, help="Comma separated list of \
            input files")
    find_range_files.add_argument("--filepath", type=str, help="Path to a root file containing \
            a list of input files")
    find_range_parser.add_argument("-loc", "--is_local", action="store_true", help='Run locally. \
            If not set will append root://cms-xrd-global.cern.ch/ \
            to the start of file names')
    find_range_parser.add_argument("--for_brilcalc", action="store_true", help='Prints the range \
            in a form compatible with the brilcalc command line tool')

    # Produce ratio config
    ratio_parser = subparsers.add_parser("produce_ratio", help="Produce ratio comparisons \
            for given numerator and denominator files")
    ratio_parser.add_argument("--numerator", type=str, required=True, help="A root file produced \
            by dijet_rdf separated by comma")
    ratio_parser.add_argument("--denominator", type=str, required=True, help="A root file \
            produced by dijet_rdf separated by comma")
    ratio_triggers = ratio_parser.add_mutually_exclusive_group()
    ratio_triggers.add_argument("--triggerlist", type=str, help="Comma separated list of triggers \
            for which plots will be produced (default value 'all')")
    ratio_triggers.add_argument("--triggerpath", type=str, help="Path to a file containing a list \
            of triggers for which plots will be produced")
    ratio_parser.add_argument("--out", type=str, required=True, default="", help="Output path \
            (output file name included)")
    ratio_parser.add_argument("--config", type=str, default="", help="Path to config file")

    # Produce responses config
    responses_parser = subparsers.add_parser("produce_responses", help="Produce responses \
            for files produced by JEC4PROMPT analysis")
    responses_files = responses_parser.add_mutually_exclusive_group(required=True)
    responses_files.add_argument("--filelist", type=str, help="Comma separated list of root files \
            produced by dijet_rdf")
    responses_files.add_argument("--filepath", type=str, help="Path to a root file containing \
            a list of output files produced by dijet_rdf")
    responses_triggers = responses_parser.add_mutually_exclusive_group()
    responses_triggers.add_argument("--triggerlist", type=str, help="Comma separated list of \
            triggers for which plots will be produced \
            (default value 'all').")
    responses_triggers.add_argument("--triggerpath", type=str, help="Path to a file containing \
            a list of triggers for which plots will be produced")
    responses_parser.add_argument("--out", type=str, default="", help="Output path")
    responses_parser.add_argument("--config", type=str, default="", help="Path to config file")

    # Produce time evolution config
    time_evolution_parser = subparsers.add_parser("produce_time_evolution", help="Produce time \
            evolution for given input files")
    time_evolution_files = time_evolution_parser.add_mutually_exclusive_group(required=True)
    time_evolution_files.add_argument("--filelist", type=str, help="Comma separated list of \
            root files produced by dijet_rdf")
    time_evolution_files.add_argument("--filepath", type=str, help="Path to a root file \
            containing a list of output files produced by dijet_rdf")
    time_evolution_triggers = time_evolution_parser.add_mutually_exclusive_group()
    time_evolution_triggers.add_argument("--triggerlist", type=str, help="Comma separated list of \
            triggers for which plots will be produced \
            (default value 'all')")
    time_evolution_triggers.add_argument("--triggerpath", type=str, help="Path to a file \
            containing a list of triggers for which plots \
            will be produced")
    time_evolution_parser.add_argument("--out", type=str, required=True, default="",
            help="Name of the output root file")
    time_evolution_parser.add_argument("--config", type=str, default="", help="Path to config file")

    # Produce vetomaps config
    vetomaps_parser = subparsers.add_parser("produce_vetomaps", help="Produce VetoMaps for files \
            produced by JEC4PROMPT analysis")
    vetomaps_files = vetomaps_parser.add_mutually_exclusive_group(required=True)
    vetomaps_files.add_argument("--filelist", type=str, help="Comma separated list of root files \
            produced by dijet_rdf")
    vetomaps_files.add_argument("--filepath", type=str, help="Path to a root file containing \
            a list of output files produced by dijet_rdf")
    vetomaps_triggers = vetomaps_parser.add_mutually_exclusive_group()
    vetomaps_triggers.add_argument("--triggerlist", type=str, help="Comma separated list of \
            triggers for which plots will be produced \
            (default value 'all')")
    vetomaps_triggers.add_argument("--triggerpath", type=str, help="Path to a file containing \
            a list of triggers for which plots will be produced")
    vetomaps_parser.add_argument("--out", type=str, default="", help="Output path")
    vetomaps_parser.add_argument("--config", type=str, default="", help="Path to config file")

    # Produce plots config
    plots_parser = subparsers.add_parser("produce_plots", help="Produce plots for given list of \
                                            input files")
    plots_files = plots_parser.add_mutually_exclusive_group(required=True)
    plots_files.add_argument("--filelist", type=str, help="Comma separated list of root files \
            produced by dijet_rdf")
    plots_files.add_argument("--filepath", type=str, help="Path to a text file containing a list \
            of output files produced by dijet_rdf")
    plots_parser.add_argument("--out", required=True, type=str, help="Output path")
    plots_parser.add_argument("--config", type=str, default="", help="Path to config file")
    plots_parser.add_argument("--all", action="store_true", 
                                help="Produce all plots in given .root files")

    # Skimming config
    skim_parser = subparsers.add_parser("skim", help="Perform skimming for\
            given list of input files")
    skim_files = skim_parser.add_mutually_exclusive_group(required=True)
    skim_files.add_argument("--filelist", type=str, help="Comma separated list of root files")
    skim_files.add_argument("--filepath", type=str, help="Path to a text file containing \
            a list of input files")
    skim_parser.add_argument("--progress_bar", action="store_true", help="Show progress bar")
    skim_parser.add_argument("--is_local", action="store_true", help='Run locally. If not set will \
            append root://cms-xrd-global.cern.ch/ to the start of file names')
    skim_triggers = skim_parser.add_mutually_exclusive_group()
    skim_triggers.add_argument('-tp', '--triggerpath', type=str, help='Path to the trigger list')
    skim_triggers.add_argument('-tl','--triggerlist', type=str, 
            help='Input files separated by commas')
    skim_parser.add_argument("--out", type=str, required=True, default="", help="Output path")
    skim_parser.add_argument("--run_range", type=str, help="Run range of the given input files \
            (run_min and run_max separated by a comma)")
    skim_parser.add_argument("--run_tag", type=str, help="Run tag")
    skim_parser.add_argument("--dataset", type=str, help="Dataset type: dijet, zjet, egamma or multijet")
    skim_parser.add_argument("--nThreads", type=int, help="Number of threads to be used \
            for multithreading")

    # Parse command line arguments, overriding config file values
    args = parser.parse_args()

    return args


#====================================================
#
# Bins for histograms and Run/IOV information
#
#====================================================
def update_run_bins(rdf: ROOT.RDF.RNode, bins: Dict) -> Dict:
    print("Recalculating run bins based on input RDF")
    run_hist = rdf.Histo1D(("run", "run", bins["runs"]["n"], bins["runs"]["bins"]), "run")
    run_hist = run_hist.GetValue()
    # Find non-zero bins
    run_bins = np.array([], dtype=float)
    last_bin = 0
    for i in range(1, run_hist.GetNbinsX()+1):
        if run_hist.GetBinContent(i) > 0:
            run_bins = np.append(run_bins, run_hist.GetBinLowEdge(i))
            last_bin = i
    run_bins = np.append(run_bins, run_hist.GetBinLowEdge(last_bin+1))
    if len(run_bins) > 1: 
        # Suspect inplace modification
        bins["runs"]["bins"] = run_bins
        bins["runs"]["n"] = len(run_bins) - 1
    return bins

def get_fill_range(IOV : str) -> tuple:
    fill_dict = {
        "Run2024H": (385814, 388000),
        "Run2024G": (383780, 385813),
        "Run2024F": (381963, 383779),
        "Run2024E": (380948, 381962),
        "Run2024D": (380253, 382947),
        "Run2024C": (379412, 380252),
        "Run2024B": (378981, 379411),
        "Run2024A": (376370, 378980),
        "Commissioning2023": (363380, 365738),
        "Run2023A": (365739, 366364),
        "Run2023B": (366365, 367079),
        "Run2023C": (367080, 369802),
        "Run2023D": (369803, 372415),
        "Run2023E": (372417, 373075),
        "Run2023F": (373076, 376371),
        "Commissioning2022": (347687, 352318),
        "Run2022A": (352319, 355064),
        "Run2022B": (355065, 355793),
        "Run2022C": (355794, 357486),
        "Run2022D": (357487, 359021),
        "Run2022E": (359022, 360331),
        "Run2022F": (360332, 362180)
    }
    if IOV not in fill_dict:
        print(f"IOV {IOV} not found in fill_dict")
        print(f"Returning default fill range (0, 1)")
        return (0, 1)

    return fill_dict[IOV]
    
def get_bins(fill_range : tuple = (376370, 380100), isMC : bool = False) -> dict:
    bins = {}
    
    # These pt and eta bins are some JEC bins
    bins["pt"] = {}
    bins["pt"]["bins"] = np.array((1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64,
        74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548,
        592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588,
        1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637,
        3832, 4037, 4252, 4477, 4713,4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000), 
        dtype=float)
    bins["pt"]["n"] = len(bins["pt"]["bins"]) - 1

    bins["eta"] = {}
    bins["eta"]["bins"] = np.array((-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839,
        -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93,
        -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957,
        -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087,
        0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218,
        1.305, 1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853,
        2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191),
        dtype=float)
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
    if not isMC:
        # bins["runs"]["bins"] = np.linspace(355065, 391370, int((391370-355065) / 10000), dtype=float)
        bins["runs"]["bins"] = np.linspace(fill_range[0], fill_range[1], 
                int((fill_range[1] - fill_range[0])), dtype=float)
        bins["runs"]["n"] = len(bins["runs"]["bins"]) - 1
    else:
        bins["runs"]["bins"] = np.linspace(0, 2, 2, dtype=float)
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
    filelist = ["Run2024C.root", "Run2024C.root", "Run2024C.root", "Run2024B.root",
            "Run2024B.root", "Run2024B.root", "Run2024A.root"]
    print(find_era(filelist))
