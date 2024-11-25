import ROOT
from processing_utils import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

from find_range import find_run_range

hist_info = [
        ("DB_direct_ratio", "Tag_pt", "DB_direct"),
        ("DB_ratio_ratio", "Tag_pt", "DB_ratio"),
        ("MPF_tag_ratio", "Tag_pt", "MPF_tag"),
        ("MPF_probe_ratio", "Probe_pt", "MPF_probe"),
        ("HDM_tag_ratio", "Tag_pt", "HDM_tag"),
        ("HDM_probe_ratio", "Probe_pt", "HDM_probe")
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
def produce_ratio(rdf_numerator, rdf_denominator, output_path, bins):
    file_ratio = ROOT.TFile.Open(f"{output_path}", "RECREATE")
    for name, x, y in hist_info:
        hn = rdf_numerator.Profile1D(("", "", bins["pt"]["n"], 
            bins["pt"]["bins"]), x, y, "weight")
        hd = rdf_denominator.Profile1D(("", "", bins["pt"]["n"],
            bins["pt"]["bins"]), x, y, "weight")
        h_ratio = hn.ProjectionX().Clone(name)
        h_ratio.Divide(hd.ProjectionX())

        h_ratio.Write()

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
        rdf_mc = (rdf_mc.Redefine("weight", f"{xsec}*genWeight"))

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    bins = get_bins()

    produce_ratio(rdf_data, rdf_mc, output_path, bins)
