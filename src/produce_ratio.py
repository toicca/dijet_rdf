import ROOT
from processing_utils import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

hist_info = [
        "multijet/MPF/MPF_multijet_PtRecoilVsResponse",
        "multijet/DB/DB_multijet_PtRecoilVsResponse"
        ]

def produce_ratio(input_file_numerator: str, input_file_denominator: str, trigger_list: List[str], output_path: str):
    file_numerator = ROOT.TFile(input_file_numerator)
    file_denominator = ROOT.TFile(input_file_denominator)
    
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers shared by the numerator and denominator files.")
        trigger_keys_numerator = {tkey.GetName() for tkey in file_numerator.GetListOfKeys()}
        trigger_keys_denominator = {tkey.GetName() for tkey in file_denominator.GetListOfKeys()}
        trigger_list = trigger_keys_numerator.intersection(trigger_keys_denominator)

    file_ratio = ROOT.TFile.Open(f"{output_path}", "RECREATE")

    for trigger_name in trigger_list:
        for namecycle in hist_info:
            hn = file_numerator.Get(f"{trigger_name}/{namecycle}")
            hd = file_denominator.Get(f"{trigger_name}/{namecycle}")
            dir_ratio = trigger_name + "/RatioResponses/" + namecycle.split("/")[0]
            if not file_ratio.GetDirectory(dir_ratio):
                file_ratio.mkdir(dir_ratio)

            if (hn.InheritsFrom("TProfile1D") and hd.InheritsFrom("TProfile1D")) \
                    or (hn.InheritsFrom("TH1D") and hd.InheritsFrom("TH1D")):
                h_name = namecycle.split("/")[-1]+"_Ratio"
                h_ratio = hn.ProjectionX().Clone(h_name)
                h_ratio.Divide(hd.ProjectionX())

                file_ratio.cd(dir_ratio)
                h_ratio.Write()
                file_ratio.cd()
            
            if hn.InheritsFrom("TH2D") and hd.InheritsFrom("TH2D"):
                h_name = namecycle.split("/")[-1]+"_Ratio"
                h_ratio = hn.Clone(h_name)
                h_ratio.Divide(hd)

                file_ratio.cd(dir_ratio)
                h_ratio.Write()
                file_ratio.cd()

            if (hn.InheritsFrom("TProfile2D") and hd.InheritsFrom("TProfile2D")) \
                    or (hn.InheritsFrom("TH3D") and hd.InheritsFrom("TH3D")):
                h_name = namecycle.split("/")[-1]+"_Ratio"
                h_ratio = hn.Clone(h_name)
                h_ratio.Divide(hd)

                file_ratio.cd(dir_ratio)
                h_ratio.Write()
                file_ratio.cd()

def run(args):
    trigger_list: List[str] = []
    files: List[str] = []
    
    if args.triggerlist:
        trigger_list = args.triggerlist.split(",")
    elif args.triggerpath:
        trigger_list = file_read_lines(args.triggerpath)

    file_numerator = args.numerator
    file_denominator = args.denominator

    output_path = args.out

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    produce_ratio(file_numerator, file_denominator, trigger_list, output_path)
