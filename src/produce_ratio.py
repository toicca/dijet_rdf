import ROOT
from processing_utils import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

from find_range import find_run_range

hist_info = [
        ("DB_direct_DataVsMC", "Tag_pt", "DB_direct"),
        ("DB_ratio_DataVsMC", "Tag_pt", "DB_ratio"),
        ("MPF_tag_DataVSMC", "Tag_pt", "MPF_tag"),
        ("MPF_probe_DataVsMC", "Probe_pt", "MPF_probe"),
        ("HDM_tag_DataVsMC", "Tag_pt", "HDM_tag"),
        ("HDM_probe_DataVsMC", "Probe_pt", "HDM_probe")
        ]

def produce_ratio(rdf_numerator, rdf_denominator, hist_config, bins):
    name = hist_config["name"]
    title = hist_config["title"]
    if hist_config["type"] == "Histo1D":
        x_bins = hist_config["x_bins"]
        x_val = hist_config["x_val"]
        hn = rdf_numerator.Histo1D((name, title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                x_val, "weight")
        hd = rdf_denominator.Histo1D((f"{name}_denom", f"{title}_denom", bins[x_bins]["n"],
            bins[x_bins]["bins"]), x_val, "weight")
        h_ratio = hn.ProjectionX().Clone(name)
        h_ratio.Divide(hd.ProjectionX())
        return h_ratio
    elif hist_config["type"] == "Profile1D":
        x_bins = hist_config["x_bins"]
        x_val = hist_config["x_val"]
        y_val = hist_config["y_val"]
        hn = rdf_numerator.Profile1D((name, title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                x_val, y_val, "weight")
        hd = rdf_denominator.Profile1D((f"{name}_denom", f"{title}_denom", bins[x_bins]["n"], bins[x_bins]["bins"]),
                x_val, y_val, "weight")
        h_ratio = hn.ProjectionX().Clone(name)
        h_ratio.Divide(hd.ProjectionX())
        return h_ratio
    else:
        raise ValueError(f"Histogram type {hist_config['type']} not supported by produce_ratio. \
                Supported types: Histo1D, Profile1D")

def run(args):
    # Shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    rdf_runs = ROOT.RDataFrame("Runs", args.data_file)
    rdf_data = ROOT.RDataFrame("Events", args.data_file)
    rdf_mc = ROOT.RDataFrame("Events", args.mc_file)
    if args.progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(rdf_runs)
        ROOT.RDF.Experimental.AddProgressBar(rdf_data)
        ROOT.RDF.Experimental.AddProgressBar(rdf_mc)

    

    min_run, max_run = find_run_range(rdf_data)

    if args.data_tag:
        output_path = f"{args.out}/J4PRatio_runs{min_run}to{max_run}_{args.data_tag}_vs_{args.mc_tag}.root"
    else:
        output_path = f"{args.out}/J4PRatio_runs{min_run}to{max_run}_vs_{args.mc_tag}.root"
   
    hist_config = read_config_file(args.hist_config)
    bins = get_bins()
   
    file_ratio = ROOT.TFile.Open(f"{output_path}", "RECREATE")
    for hist in hist_config:
        if hist.lower() == "default":
            continue
        hist = produce_ratio(rdf_data, rdf_mc, hist_config[hist], bins)
        hist.Write()
    file_ratio.Close()
