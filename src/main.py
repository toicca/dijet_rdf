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
from plotting import produce_plots, add_plots_parser
import skim

def parse_arguments():
    parser = argparse.ArgumentParser(description='JEC4PROMPT Toolkit: \
            https://github.com/toicca/dijet_rdf/tree/main')

    subparsers = parser.add_subparsers(dest='subparser_name')

    histograms.add_hist_parser(subparsers)
    find_json.add_find_json_parser(subparsers)
    find_newest.add_find_newest_parser(subparsers)
    find_range.add_find_range_parser(subparsers)
    add_plots_parser(subparsers)
    skim.add_skim_parser(subparsers)
    produce_vetomaps.add_vetomaps_parser(subparsers)
    
    # Parse command line arguments, overriding config file values
    args = parser.parse_args()

    if args.subparser_name == "skim":
        if not args.is_mc and args.mc_tag:
            raise ValueError("is_mc not set but mc_tag given")
        if args.is_mc and args.run_range:
            raise ValueError("run_range and is_mc both set")
        if (args.step is not None and args.nsteps is None) or \
                (args.nsteps is not None and args.step is None):
            raise ValueError("nsteps and step should be passed together")
        if (args.step is not None and args.nsteps is not None):
            if args.step > args.nsteps:
                raise ValueError("step should be less than nsteps")
        if args.dataset:
            print("--dataset is deprecated. Use --channel instead.")
            args.channel = args.dataset

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
