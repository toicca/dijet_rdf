import argparse
import ROOT

import find_json
import find_newest
import find_range
import histograms
import produce_ratio
import produce_responses
import produce_time_evolution
import produce_vetomaps
from plotting import produce_plots
import skim

def parse_arguments():
    parser = argparse.ArgumentParser(description='JEC4PROMPT Toolkit: \
            https://github.com/toicca/dijet_rdf/tree/main')

    subparsers = parser.add_subparsers(dest='subparser_name')

    # Histogram config: Histogram producer from skimmed files
    hist_parser = subparsers.add_parser('hist', help='Produce histograms from skimmed files.')
    hist_files = hist_parser.add_mutually_exclusive_group(required=True)
    hist_files.add_argument('-c', '--config', type=str, help='Path to the config file. If set, \
            overrides all other options.')
    hist_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
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
    find_json_parser.add_argument("--run_range", required=True, type=str, help="Run number")
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
    find_range_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
    find_range_parser.add_argument("-loc", "--is_local", action="store_true", help='Run locally. \
            If not set will append root://cms-xrd-global.cern.ch/ \
            to the start of file names')
    find_range_parser.add_argument("--for_brilcalc", action="store_true", help='Prints the range \
            in a form compatible with the brilcalc command line tool')
    find_range_parser.add_argument("--nThreads", type=int, help="Number of threads to be used \
            for multithreading")
    find_range_parser.add_argument("--progress_bar", action="store_true", help="Show progress bar")

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
    responses_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
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
    time_evolution_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
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
    vetomaps_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
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
    plots_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
    plots_parser.add_argument("--out", required=True, type=str, help="Output path")
    plots_parser.add_argument("--config", type=str, default="", help="Path to config file")
    plots_parser.add_argument("--all", action="store_true",
                                help="Produce all plots in given .root files")

    # Skimming config
    skim_parser = subparsers.add_parser("skim", help="Perform skimming for\
            given list of input files")
    skim_files = skim_parser.add_mutually_exclusive_group(required=True)
    skim_files.add_argument("--filelist", type=str, help="Comma separated list of root files")
    skim_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
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
    skim_parser.add_argument("--dataset", type=str, required=True,
            choices=["dijet", "zjet", "egamma", "multijet"],
            help="Dataset type: dijet, zjet, egamma or multijet")
    skim_parser.add_argument("--nThreads", type=int, help="Number of threads to be used \
            for multithreading")
    skim_parser.add_argument("--golden_json", type=str, help="Golden JSON for filtering")

    # Parse command line arguments, overriding config file values
    args = parser.parse_args()

    return args


if __name__ == "__main__":
    args = parse_arguments()

    command = args.subparser_name
    
    # One day when Python version >= 3.10.0
    # is used implement match here instead.
    if command == "find_json":
        find_json.run(args)
    elif command == "find_newest":
        find_newest.run(args)
    elif command == "find_range":
        find_range.run(args)
    elif command == "hist":
        histograms.run(args)
    elif command == "produce_ratio":
        produce_ratio.run(args)
    elif command ==  "produce_responses":
        produce_responses.run(args)
    elif command == "skim":
        skim.run(args)
    elif command == "produce_time_evolution":
        produce_time_evolution.run(args)
    elif command == "produce_plots":
        produce_plots.run(args)
    elif command == "produce_vetomaps":
        produce_vetomaps.run(args)
