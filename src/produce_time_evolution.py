import ROOT
from RDFHelpers import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description="Responses for dijet_rdf: https://github.com/toicca/dijet_rdf")

    files = parser.add_mutually_exclusive_group(required=True)
    files.add_argument("--filelist", type=str, help="Comma separated list of root files produced by dijet_rdf")
    files.add_argument("--filepath", type=str, help="Path to a root file containing a list of output files produced by dijet_rdf")

    triggers = parser.add_mutually_exclusive_group()
    triggers.add_argument("--triggerlist", type=str, help="Comma separated list of triggers for which plots will be produced (default value 'all')")
    triggers.add_argument("--triggerpath", type=str, help="Path to a file containing a list of triggers for which plots will be produced")

    parser.add_argument("--out", type=str, required=True, default="", help="Name of the output root file")

    parser.add_argument("--config", type=str, default="", help="Path to config file")

    args = parser.parse_args()
    
    return args

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
            xerrs = []
            for file in files:
                h = file.Get(trigger + "/" + hist)
                h.Rebin(h.GetNbinsX())
                bins.append(h.GetXaxis().GetBinCenter(1))
                vals.append(h.GetBinContent(1))
                xerrs.append(h.GetXaxis().GetBinWidth(1) / 2)
                yerrs.append(h.GetBinError(1))

            out_hist = ROOT.TGraphErrors(len(bins), np.array(bins), np.array(vals), np.array(xerrs), np.array(yerrs))
            out_hist.SetName(hist.split("/")[-1] + "_Evolution") 
            
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



if __name__ == '__main__':
    
    args = parse_arguments()
    
    trigger_list: List[str] = []
    files: List[str] = []
    
    if args.triggerpath:
        trigger_list = file_read_lines(args.triggerpath)
    elif args.triggerlist:
        trigger_list = [s.strip() for s in args.triggerlist.split(',')]

    if args.filepath:
        files = file_read_lines(args.filepath)
    else:
        print(args.filelist)
        files = [s.strip() for s in args.filelist.split(',')]

    output_path = args.out

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    # create root file, there could be an option to update
    outf = ROOT.TFile(output_path, "RECREATE")
    outf.Close()

    time_evolution(files, trigger_list, output_path) 
    cumulative_response(files, trigger_list, output_path)