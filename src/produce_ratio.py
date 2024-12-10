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

weight_info = {
    "xsec" : {
        # dijet and multijet
        "QCD-4Jets_HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 892.4,
        "QCD-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 25240000.0,
        "QCD-4Jets_HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8": 385.4,
        "QCD-4Jets_HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 126.5,
        "QCD-4Jets_HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 26.53,
        "QCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1958000.0,
        "QCD-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 96730.0,
        "QCD-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": 312200000.0,
        "QCD-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8": 13590.0,
        "QCD-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": 58840000.0,
        "QCD-4Jets_HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 3046.0,
        # zjet
        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": 6695.0,
        # egamma
        "GJ-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": 15240.0,
        "GJ-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": 8111.0,
        "GJ-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 7327.0,
        "GJ-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1541.0,
        "GJ-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 167.6,
        "GJ-4Jets_HT-600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 54.39

    },
# TODO
#    "nGenEvents" : {
#        # dijet and multijet
#        "QCD-4Jets_HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 2895970,
#        "QCD-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 5629540,
#        "QCD-4Jets_HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19537600,
#        "QCD-4Jets_HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 17527100,
#        "QCD-4Jets_HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 9212540,
#        "QCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": 18647200,
#        "QCD-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19101200,
#        "QCD-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19282700,
#        "QCD-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19122400,
#        "QCD-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1054540,
#        "QCD-4Jets_HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 18625600,
#        # zjet
#        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": 74301700
#        # egamma
#        "GJ-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-600_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???
#    }
}

def produce_ratio(rdf_numerator, rdf_denominator, hist_config, bins):
    name = hist_config["name"]
    title = hist_config["title"]
    if hist_config["type"] == "Histo1D":
        x_bins = hist_config["x_bins"]
        x_val = hist_config["x_val"]
        hn = rdf_numerator.Histo1D((name, title, bins[x_bins]["n"], bins[x_bins]["bins"]),
                x_val, "weight")
        hd = rdf_numerator.Histo1D((f"{name}_denom", f"{title}_denom", bins[x_bins]["n"],
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

    xsec = weight_info["xsec"].get(args.mc_tag)
    if xsec:
        print(f"Reweight with xsec={xsec}")
        rdf_mc = (rdf_mc.Redefine("weight", f"{xsec}*weight"))

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    
    hist_config = read_config_file(args.hist_config)
    bins = get_bins()
   
    file_ratio = ROOT.TFile.Open(f"{output_path}", "RECREATE")
    for hist in hist_config:
        if hist.lower() == "default":
            continue
        hist = produce_ratio(rdf_data, rdf_mc, hist_config[hist], bins)
        hist.Write()
    file_ratio.Close()
