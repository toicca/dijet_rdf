import ROOT
from processing_utils import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

from find_range import find_run_range
from produce_ratio import produce_ratio

def run(args):
    # Shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)
    
    files: List[str] = []
    if args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        for path in paths:
            files.extend(file_read_lines(path))
    else:
        files = [s.strip() for s in args.filelist.split(',')]

    hist_config = read_config_file(args.hist_config)

    if args.groups_of:
        if args.groups_of > len(files):
            groups = [files]
        else:
            n = args.groups_of
            groups = [files[n*i:n*i+n] for i in range(int((len(files)+n-1)/n))]
    else:
        groups = [files]
    
    events_chain = ROOT.TChain("Events")
    runs_chain = ROOT.TChain("Runs")

    for i, group in enumerate(groups):
        print(f"Processing {i+1}. group of {len(group)} files")
        events_chain = ROOT.TChain("Events")
        runs_chain = ROOT.TChain("Runs")

        for file in group:
            events_chain.Add(file)
            runs_chain.Add(file)

        rdf_data = ROOT.RDataFrame(events_chain)
        rdf_runs = ROOT.RDataFrame(runs_chain)
        rdf_mc = ROOT.RDataFrame("Events", args.mc_file)
        if args.progress_bar:
            ROOT.RDF.Experimental.AddProgressBar(rdf_runs)
            ROOT.RDF.Experimental.AddProgressBar(rdf_data)
            ROOT.RDF.Experimental.AddProgressBar(rdf_mc)

        min_run, max_run = find_run_range(rdf_data)

        if args.data_tag:
            output_path = f"{args.out}/J4PCumulative_runs{min_run}to{max_run}_{args.data_tag}_vs_{args.mc_tag}.root"
        else:
            output_path = f"{args.out}/J4PRatio_runs{min_run}to{max_run}_vs_{args.mc_tag}.root"
   
        bins = get_bins()
   
        file_ratio = ROOT.TFile.Open(f"{output_path}", "RECREATE")
        for hist in hist_config:
            if hist.lower() == "default":
                continue
            h = produce_ratio(rdf_data, rdf_mc, hist_config[hist], bins)
            h.Write()
        file_ratio.Close()
