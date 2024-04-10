import ROOT
from RDFHelpers import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

hist_names = ("PFComposition_EtaVsPhiVsProfileNEF_selected",
            "PFComposition_EtaVsPhiVsProfileCEF_selected",
            "PFComposition_EtaVsPhiVsProfileCHF_selected",
            "DB_dijet_EtaprobeVsPhiprobeVsAsymmetry",
            "Inclusive_EtaVsPhi_selected"
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

def produce_responses(file: str, trigger_list: List[str], output_path : str):
    """
    Response producer for dijet_rdf.
    """

    system = "dijet"
    methods = ["DB", "MPF"]
    bins = get_bins()
    
    file = ROOT.TFile(file, "UPDATE")
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers in the file.")
        trigger_keys = file.GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]
        
    if system == "dijet":
        pT_binLabels = ["average_Pt_dijet", "Jet_pt_tag", "Jet_pt_probe"]
    elif system == "multijet":
        pT_binLabels = ["average_Pt_multijet", "Jet_pt_lead", "pt_recoil"]
        
    for method in methods:
        for trg in trigger_list:
            path = f"{trg}/{method}/"
            if not file.GetDirectory(trg):
                file.mkdir(trg)
            if not file.GetDirectory(path):
                file.mkdir(path)
            if not file.GetDirectory(trg+"/Responses2"):
                file.mkdir(trg+"/Responses2")

            for binning in pT_binLabels:
                
                h = file.Get(path + f"{method}_{system}_PtVsEtaVsResponse_PtBin{binning}")
                # Get the projection w.r.t. y-axis
                h.GetYaxis().SetRangeUser(-2.5, 2.5)
                h2 = h.Project3D("zx")
                
                # Profile
                h3 = h2.ProfileX()
                h3.SetName(f"{method}_{system}_PtVsEtaVsResponse_PtBin{binning}_Profile")
                
                # Save
                file.cd(trg+"/Responses2")
                h2.Write()
                h3.Write()
                file.cd()
                
                if method=="MPF":
                    h = file.Get(path + f"{method}_{system}_PtVsEtaVsB_PtBin{binning}")
                else:
                    h = file.Get(path + f"{method}_{system}_PtVsEtaVsAsymmetry_PtBin{binning}")
                    
                h.GetYaxis().SetRangeUser(-2.5, 2.5)
                h2 = h.Project3D("zx")
                
                # Profile
                h3 = h2.ProfileX().ProjectionX()
                
                # Unit histogram
                unit = h3.Clone()
                unit.Reset()
                for i in range(1, unit.GetNbinsX()+1):
                    unit.SetBinContent(i, 1)
                    unit.SetBinError(i, 0)
                    
                nominator = unit.Clone()
                nominator.Add(unit, h3, c1=1.0, c2=1.0)
                denominator = unit.Clone()
                denominator.Add(unit, h3, c1=1.0, c2=-1.0)
                
                h3.Divide(nominator, denominator, 1, 1, "B")
                
                h3.SetName(f"rehti_{method}_{system}_PtVsEtaVsB_PtBin{binning}_Profile")
                
                # Save
                file.cd(trg+"/Responses2")
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
    
        
        

    