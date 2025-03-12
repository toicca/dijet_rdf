import ROOT
from utils.processing_utils import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np
import time

from find_range import find_run_range

def produce_ratio(rdf_numerator, h_denominator, hist_config, bins, i=None):
    name = hist_config["name"]
    title = hist_config["title"]
    if hist_config["type"] == "Histo1D":
        x_bins = hist_config["x_bins"]
        x_val = hist_config["x_val"]
        cut = hist_config.get("cut")
        if cut:
            hn = rdf_numerator.Filter(cut).Histo1D((f"data_{name}", title,
                bins[x_bins]["n"], bins[x_bins]["bins"]), x_val, "weight")
        else:
            hn = rdf_numerator.Histo1D((f"data_{name}", title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                    x_val, "weight")
        h_ratio = hn.Clone(f"ratio_{name}")
        h_ratio.Divide(h_denominator)
        return [hn, h_ratio]
    elif hist_config["type"] == "Profile1D":
        x_bins = hist_config["x_bins"]
        x_val = hist_config["x_val"]
        y_val = hist_config["y_val"]
        cut = hist_config.get("cut")
        if cut:
            hn = rdf_numerator.Filter(cut).Profile1D((f"data_{name}", title,
                bins[x_bins]["n"], bins[x_bins]["bins"]), x_val, y_val, "weight")
        else:
            hn = rdf_numerator.Profile1D((f"data_{name}", title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                    x_val, y_val, "weight")
        h_ratio = hn.Clone(f"ratio_{name}")
        h_ratio.Divide(h_denominator)
        return [hn, h_ratio]
    else:
        raise ValueError(f"Histogram type {hist_config['type']} not supported by produce_ratio. \
                Supported types: Histo1D, Profile1D")

def run(args):
    # Shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    print("Producing ratio(s)...")

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    mc_files = [s.strip() for s in args.mc_files.split(",")]
    chain_mc = ROOT.TChain("Events")
    for file in mc_files:
        chain_mc.Add(file)
    rdf_mc = ROOT.RDataFrame(chain_mc)

    if args.progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(rdf_mc)

    triggers = []
    if args.triggerlist:
        triggers = args.triggerlist.split(",")
    elif args.triggerpath:
        triggers = file_read_lines(args.triggerpath)

    if len(triggers) > 0:
        trg_filter = " || ".join(triggers)
        rdf_mc = (rdf_mc.Filter(trg_filter))

    data_files = [s.strip() for s in args.data_files.split(",")]
    if args.groups_of:
        if args.groups_of > len(data_files):
            groups = [data_files]
        else:
            n = args.groups_of
            groups = [data_files[n*i:n*i+n] for i in range(int((len(data_files)+n-1)/n))]
    else:
        groups = [data_files]

    bins = get_bins()
    hist_config = dict(read_config_file(args.hist_config))
    del hist_config["DEFAULT"]

    hm = {}
    for hist in hist_config:
        # Create MC histogram here early so it does not need to be
        # generated multiple times in produce_cumulative
        name = hist_config[hist]["name"]
        title = hist_config[hist]["title"]
        x_bins = hist_config[hist]["x_bins"]
        x_val = hist_config[hist]["x_val"]
        y_val = hist_config[hist]["y_val"]
        cut = hist_config[hist].get("cut")
        if cut:
            h = rdf_mc.Filter(cut).Profile1D((f"mc_{name}", title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                x_val, y_val, "weight")
        else:
            h = rdf_mc.Profile1D((f"mc_{name}", title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                x_val, y_val, "weight")
        hm[hist] = h.ProjectionX()

    lm = {}
    for hist in hist_config:
        x_val = hist_config[hist]["x_val"]
        y_val = hist_config[hist]["y_val"]
        cut = hist_config[hist].get("cut")
        if cut:
            lm[hist] = rdf_mc.Filter(cut).Stats(y_val, "weight")
        else:
            lm[hist] = rdf_mc.Stats(y_val, "weight")

    for i, group in enumerate(groups):
        chain_data = ROOT.TChain("Events")
        chain_runs = ROOT.TChain("Runs")

        lds = []
        for j, file in enumerate(group):
            chain_data.Add(file)
            chain_runs.Add(file)
            rdf = ROOT.RDataFrame("Events", file)

        rdf_data = ROOT.RDataFrame(chain_data)
        rdf_runs = ROOT.RDataFrame(chain_runs)

        if args.progress_bar:
            ROOT.RDF.Experimental.AddProgressBar(rdf_runs)
            ROOT.RDF.Experimental.AddProgressBar(rdf_data)

        if len(triggers) > 0:
            trg_filter = " || ".join(triggers)
            rdf_data = (rdf_data.Filter(trg_filter))

        min_run = int(rdf_data.Min("min_run").GetValue())
        max_run = int(rdf_data.Max("max_run").GetValue())
        print(f"Group run range: [{min_run}, {max_run}]")
        if args.data_tag:
            output_path = f"{args.out}/J4PRatio_runs{min_run}to{max_run}_{args.data_tag}_vs_{args.mc_tag}.root"
        else:
            output_path = f"{args.out}/J4PRatio_runs{min_run}to{max_run}_vs_{args.mc_tag}.root"

        hists_out = []
        for hist in hist_config:
            [hn, h] = produce_ratio(rdf_data, hm[hist], hist_config[hist], bins)
            hists_out.append(hm[hist])
            hists_out.append(hn)
            hists_out.append(h)
            print(f"{hist} processed")

        file_ratio = ROOT.TFile.Open(f"{output_path}", "RECREATE")
        for h in hists_out:
            h.Write()
        file_ratio.Close()

        if len(groups) > 1:
            print(f"{i+1}/{len(groups)} groups processed")
