import ROOT
from processing_utils import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

def time_evolution(filelist: List[str], trigger_list: List[str], output_path: str):
    
    time_histograms = ("multijet/MPF/MPF_multijet_RunVsResponse",
                       "multijet/DB/DB_multijet_RunVsResponse",)

    # Load the root files to be analyzed
    files = [ROOT.TFile(file, "READ") for file in filelist]

    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers from the first file.")
        trigger_keys = files[0].GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]

    out_file = ROOT.TFile(output_path, "UPDATE")
    for trigger in trigger_list:
        for hist in time_histograms:
            bins = []
            vals = []
            yerrs = []
            upper_edges = []
            for file in files:
                h = file.Get(trigger + "/" + hist)
                h.Rebin(h.GetNbinsX())
                bins.append(h.GetBinLowEdge(1))
                vals.append(h.GetBinContent(1))
                yerrs.append(h.GetBinError(1))
                upper_edges.append(h.GetBinWidth(1) + h.GetBinLowEdge(1))

            # Sort the bins and values by the bin edges
            bins, vals, yerrs, upper_edges = zip(*sorted(zip(bins, vals, yerrs, upper_edges)))
            bins = list(bins) + [upper_edges[-1]]

            out_hist = ROOT.TH1D(hist.split("/")[-1] + "_Evolution", hist.split("/")[-1] + "_Evolution", len(bins)-1, np.array(bins, dtype=float))

            for i in range(len(vals)):
                out_hist.SetBinContent(i+1, vals[i])
                out_hist.SetBinError(i+1, yerrs[i])
            
            if not out_file.GetDirectory(trigger + "/TimeEvolution/" + hist.split("/")[0]):
                out_file.mkdir(trigger+"/TimeEvolution/"+hist.split("/")[0])
            out_file.cd(trigger+"/TimeEvolution/"+hist.split("/")[0])
            out_hist.Write()
            out_file.cd()

    out_file.Close()


def cumulative_response(filelist: List[str], trigger_list: List[str], output_path: str):
    
    response_histograms = ("multijet/MPF/MPF_multijet_PtRecoilVsResponseCorrected",
                           "multijet/DB/DB_multijet_PtRecoilVsResponseCorrected",)
    
    # Load the root files to be analyzed
    files = [ROOT.TFile(file, "READ") for file in filelist]

    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers from the first file.")
        trigger_keys = files[0].GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]

    out_file = ROOT.TFile(output_path, "UPDATE")

    bins = get_bins()

    for trigger in trigger_list:
        for hist in response_histograms:
            cum_prof = ROOT.TProfile(hist.split("/")[-1] + "_Cumulative", hist.split("/")[-1] + "_Cumulative", bins["pt"]["n"], bins["pt"]["bins"])

            for file in files:
                h = file.Get(trigger + "/" + hist)
                cum_prof.Add(h)

            if not out_file.GetDirectory(trigger + "/CumulativeResponses/" + hist.split("/")[0]):
                out_file.mkdir(trigger + "/CumulativeResponses/" + hist.split("/")[0])
            out_file.cd(trigger + "/CumulativeResponses/" + hist.split("/")[0])
            cum_prof.Write()
            out_file.cd()



def run(args):
    trigger_list: List[str] = []
    files: List[str] = []
    
    if args.triggerlist:
        trigger_list = args.triggerlist.split(",")
    elif args.triggerpath:
        trigger_list = file_read_lines(args.triggerpath)

    # Split the file list and trigger list if they are given as a string
    if args.filelist:
        files= args.filelist.split(",")
    elif args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        for path in paths:
            files.extend(file_read_lines(path, find_ROOT=True))
    else:
        raise ValueError("No file list provided")

    output_path = args.out

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    # create root file, there could be an option to update
    outf = ROOT.TFile(output_path, "RECREATE")
    outf.Close()

    time_evolution(files, trigger_list, output_path) 
    cumulative_response(files, trigger_list, output_path)
