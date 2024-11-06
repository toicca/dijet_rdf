import ROOT
import configparser
import os
import subprocess
import json
# import tomllib
from processing_utils import find_site, get_bins
from RDFHelpers import read_config_file

def create_histogram(rdf, hist_config, bins):
    if hist_config["type"] == "Histo1D":
        return rdf.Histo1D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"]),
                            hist_config["x_val"], "weight")
    elif hist_config["type"] == "Histo2D":
        return rdf.Histo2D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], "weight")
    elif hist_config["type"] == "Histo3D":
        return rdf.Histo3D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"],
                            bins[hist_config["z_bins"]]["n"], bins[hist_config["z_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], hist_config["z_val"], "weight")
    elif hist_config["type"] == "Profile1D":
        return rdf.Profile1D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], "weight")
    elif hist_config["type"] == "Profile2D":
        return rdf.Profile2D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], hist_config["z_val"], "weight")
    elif hist_config["type"] == "Profile3D":
        return rdf.Profile3D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"],
                            bins[hist_config["z_bins"]]["n"], bins[hist_config["z_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], hist_config["z_val"], "weight")
    else:
        raise ValueError(f"Unknown histogram type: {hist_config['type']}")

def make_histograms(args):
    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    events_chain = ROOT.TChain("Events")
    runs_chain = ROOT.TChain("Runs")
                
    # Split the file list and trigger list if they are given as a string
    if args.filelist:
        filelist = args.filelist.split(",")
    elif args.filepath:
        filelist = file_read_lines(args.filepath, find_ROOT=True)
    else:
        raise ValueError("No file list provided")

    # Load the files
    for file in filelist:
        if not args.is_local:
            events_chain.Add(f"root://cms-xrd-global.cern.ch/{file}")
            runs_chain.Add(f"root://cms-xrd-global.cern.ch/{file}")
        else:
            events_chain.Add(file)
            runs_chain.Add(file)

    events_rdf = ROOT.RDataFrame(events_chain)
    runs_rdf = ROOT.RDataFrame(runs_chain)

    if args.progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(events_rdf)

    # with open(config['histogram_config'], 'rb') as f:
        # hist_config = tomllib.load(f)
    hist_config = configparser.ConfigParser()
    hist_config.read(args.hist_config)
    hist_config = dict(hist_config)
    bins = get_bins()

    histograms = {}

    for hist in hist_config:
        if hist.lower() == "default":
            continue
        histograms[hist] = create_histogram(events_rdf, hist_config[hist], bins).GetValue()

    return histograms

def get_values(histograms):
    values = {}
    for hist in histograms:
        values[hist] = histograms[hist].GetValue()
    return values

def save_histograms(histograms, args):
    run_range = args.run_range.split(",")
    assert(len(run_range) == 2)

    if args.out:
        if not os.path.exists(args.out):
            os.makedirs(args.out)
        output_file = ROOT.TFile(
                f"{args.out}/J4PHists_runs{run_range[0]}to{run_range[1]}_{args.run_tag}.root", 
                "RECREATE")
    else:
        output_file = ROOT.TFile(
                f"J4PHists_runs{run_range[0]}to{run_range[1]}_{args.run_tag}.root",
                "RECREATE")

    for hist in histograms:
        histograms[hist].Write()

    output_file.Close()

def run(args):
    """
    config = {
        'filelist': ['J4PSkim_runs379413to379415_20240924.root'],
        'is_local': True,
        'progress_bar': True,
        'histogram_config': 'histograms.toml',
        'run_range': (379413, 379415),
        'run_tag': '20240920',
        'nThreads': 8
    }
    """

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

    histograms = make_histograms(args)
    save_histograms(histograms, args)    
