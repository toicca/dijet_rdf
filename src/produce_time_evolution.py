import ROOT
from utils.processing_utils import file_read_lines, read_config_file, get_bins
from typing import List
import json
import argparse, configparser
import numpy as np
import time

from find_range import find_run_range

def update_state(state):
    add_produce_time_evolution_parser(state)
    state.valfuncs["produce_time_evolution"] = validate_args
    state.commands["produce_time_evolution"] = run

def add_produce_time_evolution_parser(state):
    subparsers = state.subparsers
    time_evolution_parser = subparsers.add_parser("produce_time_evolution", help="Produce time \
            evolution for given input files")
    time_evolution_files = time_evolution_parser.add_mutually_exclusive_group(required=True)
    time_evolution_files.add_argument("--filelist", type=str, help="Comma separated list of \
            input files")
    time_evolution_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
    time_evolution_parser.add_argument("-tf", "--triggerfile", type=str, help="Path to the .json file \
        containing the list of triggers to be used.")
    time_evolution_parser.add_argument("-ch", "--channel", type=str, help="Channel to be used",
                                       choices=state.channels)
    time_evolution_parser.add_argument("--out", type=str, required=True, default="", help="Output path \
            (output file name included)")
    time_evolution_parser.add_argument("--data_tag", type=str, help="data tag")
    time_evolution_parser.add_argument("--nThreads", type=int, help="Number of threads to be used \
            for multithreading")
    time_evolution_parser.add_argument("--progress_bar", action="store_true", help="Show progress bar")
    time_evolution_parser.add_argument('-hconf', '--hist_config', required=True, type=str, help='Path to the histogram \
            config file.')

def validate_args(args):
    pass

def data_hists(rdf, hist_config, bins):
    hd = {}
    hd["int_lumi"] = rdf.Mean("int_lumi").GetValue()
    hd["min_run"] = rdf.Min("min_run").GetValue()
    hd["max_run"] = rdf.Max("max_run").GetValue()
    for hist in hist_config:
        name = hist_config[hist]["name"]
        title = hist_config[hist]["title"]
        if hist_config[hist]["type"] == "Histo1D":
            x_bins = hist_config[hist]["x_bins"]
            x_val = hist_config[hist]["x_val"]
            cut = hist_config[hist].get("cut")
            if cut:
                hd[hist] = rdf.Filter(cut).Histo1D((name, title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                        x_val, "weight")
            else:
                hd[hist] = rdf.Histo1D((name, title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                    x_val, "weight")
        elif hist_config[hist]["type"] == "Profile1D":
            x_bins = hist_config[hist]["x_bins"]
            x_val = hist_config[hist]["x_val"]
            y_val = hist_config[hist]["y_val"]
            cut = hist_config[hist].get("cut")
            if cut:
                hd[hist] = rdf.Filter(cut).Profile1D((name, title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                        x_val, y_val, "weight")
            else:
                hd[hist] = rdf.Profile1D((name, title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                        x_val, y_val, "weight")
    return hd

def lumi_data(rdf, hist_config, triggers):
    ld = {}
    ld["int_lumi"] = rdf.Mean("int_lumi")
    ld["min_run"] = rdf.Min("min_run")
    ld["max_run"] = rdf.Max("max_run")

    trg_filter = "1"
    if len(triggers) > 0:
        trg_filter = " || ".join(triggers)

    for hist in hist_config:
        x_val = hist_config[hist]["x_val"]
        y_val = hist_config[hist]["y_val"]
        cut = hist_config[hist].get("cut")
        if cut:
            ld[hist] = rdf.Filter(trg_filter).Filter(cut).Stats(y_val, "weight")
        else:
            ld[hist] = rdf.Filter(trg_filter).Stats(y_val, "weight")
    return ld


def produce_time_evolution(lds, hist_config, bins, logging):
    handles = []
    for ld in lds:
        for hist in hist_config:
            handles.append(ld[hist])
    ROOT.RDF.RunGraphs(handles)

    clumi = 0.0
    lumi = 0.0
    lumi_bins = [0.0]
    for ld in lds:
        lumi = ld["int_lumi"].GetValue()
        clumi += lumi
        lumi_bins.append(clumi)

    hs = {}
    hs["min_runs"] = ROOT.TH1D(f"min_runs", "min_runs", len(lumi_bins)-1, np.array(lumi_bins))
    hs["max_runs"] = ROOT.TH1D(f"max_runs", "max_runs", len(lumi_bins)-1, np.array(lumi_bins))
    for hist in hist_config:
        name = hist_config[hist]["name"]
        title = f"{hist};Cumulative luminosity (fb^{{-1}});<Response>"

        hs[hist] = ROOT.TH1D(f"{name}_int_lumi", title, len(lumi_bins)-1,
                np.array(lumi_bins))

    logging.info("Filling time evolution histograms")
    start_time_evolution = time.time()
    for i, ld in enumerate(lds):
        int_lumi = ld["int_lumi"].GetValue()
        min_run = ld["min_run"].GetValue()
        max_run = ld["max_run"].GetValue()

        hs["min_runs"].SetBinContent(i+1, min_run)
        hs["max_runs"].SetBinContent(i+1, max_run)
        for hist in hist_config:
            start_ratio = time.time()
            md = ld[hist].GetValue().GetMean()
            ed = ld[hist].GetValue().GetMeanErr()
            hs[hist].SetBinContent(i+1, md)
            hs[hist].SetBinError(i+1, ed)
    logging.info(f"Done filling time evolution histograms (took {time.time() - start_time_evolution} s)")
    hists = [hs["min_runs"], hs["max_runs"]]
    for hist in hist_config:
        hists.append(hs[hist])

    return hists

def run(state):
    args = state.args
    logging = state.logger
    # Shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    logging.info("Producing time evolution histograms...")

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    files = []
    if args.filelist:
        files = [s.strip() for s in args.filelist.split(",")]
    elif args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        for path in paths:
            files.extend(file_read_lines(path, find_ROOT=True))
    else:
        raise ValueError("No file list provided")

    with open(args.triggerfile, "r") as f:
        triggers = json.load(f)[args.channel]

    for trigger in triggers:
        triggers[trigger] = triggers[trigger]["cut"]

    bins = get_bins()
    hist_config = dict(read_config_file(args.hist_config))
    del hist_config["DEFAULT"]

    chain_data = ROOT.TChain("Events")
    chain_runs = ROOT.TChain("Runs")

    lds = []
    for file in files:
        rdf = ROOT.RDataFrame("Events", file)

        ld = lumi_data(rdf, hist_config, triggers)
        lds.append(ld)

    min_run = int(min([ld["min_run"].GetValue() for ld in lds]))
    max_run = int(max([ld["max_run"].GetValue() for ld in lds]))
    logging.info(f"Group run range: [{min_run}, {max_run}]")
    if args.data_tag:
        output_path = f"{args.out}/J4PTimeEvolution_runs{min_run}to{max_run}_{args.data_tag}.root"
    else:
        output_path = f"{args.out}/J4PTimeEvolution_runs{min_run}to{max_run}.root"

    output_file = ROOT.TFile.Open(f"{output_path}", "RECREATE")
    hs = produce_time_evolution(lds, hist_config, bins, logging)
    for h in hs:
            h.Write()
    output_file.Close()
