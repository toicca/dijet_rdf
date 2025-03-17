import argparse
import configparser
import json
import os

import ROOT

# import tomllib
from utils.processing_utils import file_read_lines, get_bins, read_config_file


def update_state(state):
    add_hist_parser(state.subparsers)
    state.valfuncs["hist"] = validate_args
    state.commands["hist"] = run


def add_hist_parser(subparsers):
    hist_parser = subparsers.add_parser(
        "hist", help="Produce histograms from skimmed files."
    )

    # Input files
    hist_files = hist_parser.add_mutually_exclusive_group(required=True)
    hist_files.add_argument(
        "-c",
        "--config",
        type=str,
        help="Path to the config file. If set, \
            overrides all other options.",
    )
    hist_files.add_argument(
        "-fp",
        "--filepaths",
        type=str,
        help="Comma separated list of \
            text files containing input files (one input file per line).",
    )
    hist_files.add_argument(
        "-fl", "--filelist", type=str, help="Input files separated by commas."
    )

    # Other
    hist_parser.add_argument(
        "-tf",
        "--triggerfile",
        type=str,
        help="Path to the .json file containing \
            triggers.",
    )
    hist_parser.add_argument(
        "-ch", "--channel", type=str, help="Channel associated with the triggers."
    )
    hist_parser.add_argument(
        "-reg",
        "--regions",
        type=str,
        help="Comma separated list of \
            .ini files with cuts for different regions.",
    )
    hist_parser.add_argument(
        "-loc",
        "--is_local",
        action="store_true",
        help="Run locally. If not \
            set will append root://cms-xrd-global.cern.ch/ \
            to the start of file names.",
    )
    hist_parser.add_argument(
        "-pbar", "--progress_bar", action="store_true", help="Show progress bar."
    )
    hist_parser.add_argument(
        "-hconf",
        "--hist_config",
        type=str,
        help="Path to the histogram \
            config file.",
    )
    hist_parser.add_argument(
        "--run_range",
        type=str,
        help="Run range of the given input files \
            (run_min and run_max separated by a comma)",
    )
    hist_parser.add_argument("--run_tag", type=str, help="Run tag")
    hist_parser.add_argument(
        "--nThreads",
        type=int,
        help="Number of threads to be used \
            for multithreading",
    )
    hist_parser.add_argument(
        "--out",
        type=str,
        required=True,
        default="",
        help="Output path \
            (output file name included)",
    )


def validate_args(args):
    if args.triggerfile and not args.channel:
        raise argparse.ArgumentError(
            None, "--channel argument is required when --triggerfile is provided"
        )


def create_histogram(rdf, hist_config, bins, triggers):
    if len(triggers) > 0:
        trg_filter = ""
        for trigger in triggers:
            trg_filter += f"({triggers[trigger]}) || "

        trg_filter = trg_filter[:-4]
        rdf = rdf.Filter(trg_filter)

    cut = hist_config.get("cut")

    if cut:
        rdf = rdf.Filter(cut)

    if hist_config["type"] == "Histo1D":
        return rdf.Histo1D(
            (
                hist_config["name"],
                hist_config["title"],
                bins[hist_config["x_bins"]]["n"],
                bins[hist_config["x_bins"]]["bins"],
            ),
            hist_config["x_val"],
            "weight",
        )
    elif hist_config["type"] == "Histo2D":
        return rdf.Histo2D(
            (
                hist_config["name"],
                hist_config["title"],
                bins[hist_config["x_bins"]]["n"],
                bins[hist_config["x_bins"]]["bins"],
                bins[hist_config["y_bins"]]["n"],
                bins[hist_config["y_bins"]]["bins"],
            ),
            hist_config["x_val"],
            hist_config["y_val"],
            "weight",
        )
    elif hist_config["type"] == "Histo3D":
        return rdf.Histo3D(
            (
                hist_config["name"],
                hist_config["title"],
                bins[hist_config["x_bins"]]["n"],
                bins[hist_config["x_bins"]]["bins"],
                bins[hist_config["y_bins"]]["n"],
                bins[hist_config["y_bins"]]["bins"],
                bins[hist_config["z_bins"]]["n"],
                bins[hist_config["z_bins"]]["bins"],
            ),
            hist_config["x_val"],
            hist_config["y_val"],
            hist_config["z_val"],
            "weight",
        )
    elif hist_config["type"] == "Profile1D":
        return rdf.Profile1D(
            (
                hist_config["name"],
                hist_config["title"],
                bins[hist_config["x_bins"]]["n"],
                bins[hist_config["x_bins"]]["bins"],
            ),
            hist_config["x_val"],
            hist_config["y_val"],
            "weight",
        )
    elif hist_config["type"] == "Profile2D":
        return rdf.Profile2D(
            (
                hist_config["name"],
                hist_config["title"],
                bins[hist_config["x_bins"]]["n"],
                bins[hist_config["x_bins"]]["bins"],
                bins[hist_config["y_bins"]]["n"],
                bins[hist_config["y_bins"]]["bins"],
            ),
            hist_config["x_val"],
            hist_config["y_val"],
            hist_config["z_val"],
            "weight",
        )
    elif hist_config["type"] == "Profile3D":
        return rdf.Profile3D(
            (
                hist_config["name"],
                hist_config["title"],
                bins[hist_config["x_bins"]]["n"],
                bins[hist_config["x_bins"]]["bins"],
                bins[hist_config["y_bins"]]["n"],
                bins[hist_config["y_bins"]]["bins"],
                bins[hist_config["z_bins"]]["n"],
                bins[hist_config["z_bins"]]["bins"],
            ),
            hist_config["x_val"],
            hist_config["y_val"],
            hist_config["z_val"],
            "weight",
        )
    else:
        raise ValueError(f"Unknown histogram type: {hist_config['type']}")


def make_histograms(args, logger):
    bins = get_bins()

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    events_chain = ROOT.TChain("Events")
    runs_chain = ROOT.TChain("Runs")

    # Split the file list and trigger list if they are given as a string
    if args.filelist:
        filelist = [s.strip() for s in args.filelist.split(",")]
    elif args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        filelist = []
        for path in paths:
            filelist.extend(file_read_lines(path, find_ROOT=True))
    else:
        raise ValueError("No file list provided")

    # Load the trigger json
    with open(args.triggerfile) as f:
        triggers = json.load(f)[args.channel]

    for trigger in triggers:
        triggers[trigger] = triggers[trigger]["cut"]

    # Load the files
    for file in filelist:
        if not args.is_local:
            events_chain.Add(f"{args.redirector}{file}")
            runs_chain.Add(f"{args.redirector}{file}")
        else:
            events_chain.Add(file)
            runs_chain.Add(file)

    events_rdf = ROOT.RDataFrame(events_chain)

    if args.progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(events_rdf)

    region_configs = {}
    if args.regions:
        regions = args.regions.split(",")
        for region in regions:
            region_config = configparser.ConfigParser()
            region_config.read(region)
            region_config = dict(region_config)
            for region in region_config:
                region_configs[region] = region_config[region].get("cut")

    hist_configs = args.hist_config.split(",")
    histograms = {}
    all_hists = []
    for hist_in in hist_configs:
        hist_config = configparser.ConfigParser()
        hist_config.read(hist_in)
        hist_config = dict(hist_config)

        for hist in hist_config:
            for region in region_configs:
                hist_config[hist + "_" + region] = hist_config[hist].copy()
                hist_config[hist + "_" + region]["cut"] = (
                    hist_config[hist + "_" + region]["cut"]
                    + " && "
                    + region_configs[region]
                )

        for hist in hist_config:
            if hist.lower() == "default":
                continue
            all_hists.append(
                create_histogram(events_rdf, hist_config[hist], bins, triggers)
            )

    for hist in all_hists:
        histograms[hist.GetName()] = hist.GetValue()

    return histograms


def get_values(histograms):
    values = {}
    for hist in histograms:
        values[hist] = histograms[hist].GetValue()
    return values


def save_histograms(histograms, args, logger):
    range_str = ""
    if args.run_range:
        run_range = args.run_range.split(",")
        assert len(run_range) == 2
        range_str = f"runs{run_range[0]}to{run_range[1]}_"

    if args.out:
        if not os.path.exists(args.out):
            os.makedirs(args.out)
        output_file = ROOT.TFile(
            f"{args.out}/J4PHists_{range_str}{args.run_tag}.root", "RECREATE"
        )
    else:
        output_file = ROOT.TFile(f"J4PHists_{range_str}{args.run_tag}.root", "RECREATE")

    for hist in histograms:
        histograms[hist].Write()

    output_file.Close()
    logger.info(f"Histograms saved to {output_file.GetName()}")


def run(state):
    args = state.args
    logger = state.logger

    # shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    # If config file is set, override all arguments with ones set in there
    if args.config:
        config = read_config_file(args.config)
        for arg, value in config["GENERAL"].items():
            # Do a type conversion for the option
            if arg == "is_local" or arg == "nThreads" or arg == "progress_bar":
                if value == "":
                    setattr(args, arg, 0)
                else:
                    setattr(args, arg, int(value))
            else:
                setattr(args, arg, value)

    histograms = make_histograms(args, logger)
    save_histograms(histograms, args, logger)
