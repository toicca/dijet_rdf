import ROOT
from RDFHelpers import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

response_histos = (("multijet", "MPF", "MPF_multijet_PtAvgVsEtaVsResponse"),
                    ("multijet", "MPF", "MPF_multijet_PtRecoilVsEtaVsResponse"),
                    ("multijet", "MPF", "MPF_multijet_PtLeadVsEtaVsResponse"),
                    ("multijet", "DB", "DB_multijet_PtAvgVsEtaVsResponse"),
                    ("multijet", "DB", "DB_multijet_PtRecoilVsEtaVsResponse"),
                    ("multijet", "DB", "DB_multijet_PtLeadVsEtaVsResponse"),
                    ("dijet", "MPF", "MPF_dijet_PtAvgVsEtaVsResponse"),
                    ("dijet", "MPF", "MPF_dijet_PtProbeVsEtaVsResponse"),
                    ("dijet", "MPF", "MPF_dijet_PtTagVsEtaVsResponse"),
                    ("dijet", "DB", "DB_dijet_PtAvgVsEtaVsResponse"),
                    ("dijet", "DB", "DB_dijet_PtProbeVsEtaVsResponse"),
                    ("dijet", "DB", "DB_dijet_PtTagVsEtaVsResponse"),
)
resolution_histos = (("dijet", "DB", "DB_dijet_PtAvgVsEtaVsA"),
                    ("dijet", "DB", "DB_dijet_PtProbeVsEtaVsA"),
                    ("dijet", "DB", "DB_dijet_PtTagVsEtaVsA"),
)

derived_histos = (("multijet", "MPF", "MPF_multijet_PtAvgVsEtaVsB"),
                    ("multijet", "MPF", "MPF_multijet_PtRecoilVsEtaVsB"),
                    ("multijet", "MPF", "MPF_multijet_PtLeadVsEtaVsB"),
                    ("multijet", "DB", "DB_multijet_PtAvgVsEtaVsA"),
                    ("multijet", "DB", "DB_multijet_PtRecoilVsEtaVsA"),
                    ("multijet", "DB", "DB_multijet_PtLeadVsEtaVsA"),
)

def parse_arguments():
    parser = argparse.ArgumentParser(description="Dijet responses for dijet_rdf: https://github.com/toicca/dijet_rdf")

    files = parser.add_mutually_exclusive_group(required=True)
    files.add_argument("--filelist", type=str, help="Comma separated list of root files produced by dijet_rdf")
    files.add_argument("--filepath", type=str, help="Path to a root file containing a list of output files produced by dijet_rdf")

    triggers = parser.add_mutually_exclusive_group()
    triggers.add_argument("--triggerlist", type=str, help="Comma separated list of triggers for which plots will be produced (default value 'all')")
    triggers.add_argument("--triggerpath", type=str, help="Path to a file containing a list of triggers for which plots will be produced")

    parser.add_argument("--out", type=str, default="", help="Output path")

    parser.add_argument("--config", type=str, default="", help="Path to config file")

    args = parser.parse_args()
    
    return args

def produce_resolutions(file: str, trigger_list: List[str], output_path : str):
    """
    Resolution producer for dijet_rdf.
    """

    bins = get_bins()
    
    file = ROOT.TFile(file, "UPDATE")
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers in the file.")
        trigger_keys = file.GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]
        
    for trg in trigger_list:
        for system, method, histogram in resolution_histos:
            path = f"{trg}/{system}/{method}/"
            resolution_path = f"{trg}/{system}/Resolutions"
            if not file.GetDirectory(path):
                file.mkdir(path)
            if not file.GetDirectory(resolution_path):
                file.mkdir(resolution_path)

            h = file.Get(path + histogram)
            h.RebinY(4)
            resolutions = ROOT.TH2D("resolutions_"+histogram.replace("VsA", "").replace("VsEta", "").replace(f"_{method}_{system}", ""), h.GetTitle(),
                                    bins["pt"]["n"], bins["pt"]["bins"], h.GetNbinsY(), h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax())
            projections = ROOT.TH2D("projections_"+histogram.replace("VsA", "").replace("VsEta", "").replace(f"_{method}_{system}", ""), h.GetTitle(),
                                    bins["pt"]["n"], bins["pt"]["bins"], bins["pt"]["n"], bins["pt"]["bins"])


            # TH3.FitSlicesZ?
            for i in range(1, h.GetNbinsX()+1):
                for j in range(1, h.GetNbinsY()+1):
                    # Fit to the asymmetry in the bin
                    proj = h.ProjectionZ(f"proj_{i}_{j}", i, i, j, j)
                    # Find the mean
                    mean = proj.GetMean()
                    # Find the sigma
                    sigma = 0
                    for k in range(1, proj.GetNbinsX()+1):
                        sigma += (proj.GetBinContent(k) - mean)**2
                    sigma = proj.GetRMS()

                    if sigma > 0: 
                        gaussian = ROOT.TF1("my_gaus", "gaus", mean - sigma, mean + sigma)
                        proj.Fit(gaussian, "Q L N")
                        fsigma = gaussian.GetParameter(2)
                        resolutions.SetBinContent(i, j, fsigma)
                        resolutions.SetBinError(i, j, gaussian.GetParError(2))
                    else:
                        resolutions.SetBinContent(i, j, 0)
                        resolutions.SetBinError(i, j, 0)
                    # Save the projection
                    for k in range(1, proj.GetNbinsX()+1):
                        projections.SetBinContent(i, k, proj.GetBinContent(k))

            resolutions.SetName(method + "_projected_resolution_"+histogram.replace("VsA", "").replace("VsEta", "").replace(f"_{method}_{system}", ""))
            
            # Save
            file.cd(resolution_path)
            # h2.Write()
            resolutions.Write()
            file.cd()


def produce_responses(file: str, trigger_list: List[str], output_path : str):
    """
    Response producer for dijet_rdf.
    """

    bins = get_bins()
    
    file = ROOT.TFile(file, "UPDATE")
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers in the file.")
        trigger_keys = file.GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]
        
    for trg in trigger_list:
        for system, method, histogram in response_histos:
            path = f"{trg}/{system}/{method}/"
            response_path = f"{trg}/{system}/Responses"
            if not file.GetDirectory(path):
                file.mkdir(path)
            if not file.GetDirectory(response_path):
                file.mkdir(response_path)

            h = file.Get(path + histogram)
            # Get the projection w.r.t. y-axis
            h.GetYaxis().SetRangeUser(-2.5, 2.5)
            h2 = h.Project3D("zx")
            
            # Profile
            h3 = h2.ProfileX()
            h3.SetName(method + "_projected_response_"+histogram.replace("VsResponse", "").replace("VsEta", "").replace(f"_{method}_{system}", ""))
            
            # Save
            file.cd(response_path)
            # h2.Write()
            h3.Write()
            file.cd()
            
        for system, method, histogram in derived_histos:    
            path = f"{trg}/{system}/{method}/"
            response_path = f"{trg}/{system}/Responses"
            if not file.GetDirectory(path):
                file.mkdir(path)
            if not file.GetDirectory(response_path):
                file.mkdir(response_path)

            h = file.Get(path + histogram)
            # Get the projection w.r.t. y-axis
            h.GetYaxis().SetRangeUser(-2.5, 2.5)
            h2 = h.Project3D("zx")
            
            # # Profile
            h3 = h2.ProfileX().ProjectionX()
            
            # # Unit histogram
            unit = h3.Clone()
            unit.Reset()
            for i in range(1, unit.GetNbinsX()+1):
                unit.SetBinContent(i, 1)
                unit.SetBinError(i, 0)
                
            nominator = unit.Clone()
            nominator.Add(unit, h3, c1=1.0, c2=1.0)
            denominator = unit.Clone()
            denominator.Add(unit, h3, c1=1.0, c2=-1.0)
            
            h3.Divide(nominator, denominator, 1, 1)
            
            h3.SetName(method + "_derived_response_"+histogram.replace("VsResponse", "").replace("VsEta", "").replace(f"_{method}_{system}", ""))
            
            # Save
            file.cd(response_path)
            h3.Write()
        


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
        files = [s.strip() for s in args.filelist.split(',')]

    output_path = args.out

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    for file in files:
        produce_responses(file, trigger_list, output_path)
        produce_resolutions(file, trigger_list, output_path)
        
        

    