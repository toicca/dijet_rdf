import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import configparser
from typing import List
import pathlib
import cmsstyle as CMS
import time

def file_read_lines(file: str) -> List[str]:
    with open(file) as f:
        return [line.strip() for line in f.readlines()]

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot producer for dijet_rdf: https://github.com/toicca/dijet_rdf")

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

def read_config_file(file: str) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.read(file)

    return config


def produce_plots(file, output_path, config, trigger_list=[]):
    triggers = trigger_list
    if len(triggers) == 0:
        trigger_keys = file.GetListOfKeys()
        triggers = [tkey.GetName() for tkey in trigger_keys]

    for trigger in triggers:
        systems = file.Get(trigger)
        for system_key in systems.GetListOfKeys():
            system_name = system_key.GetName()
            methods = systems.Get(system_name)
            for method_key in methods.GetListOfKeys():
                method_name = method_key.GetName()
                pathlib.Path("{}/{}/{}/{}".format(output_path, trigger, system_name, method_name)).mkdir(exist_ok=True, parents=True)
                hists = methods.Get(method_name)
                for hist_key in hists.GetListOfKeys():
                    hist_name = hist_key.GetName()
                    hist = hists.Get(hist_name)

                    # ignore 3D histograms and profiles for now
                    if hist.InheritsFrom("TH3D") or hist.InheritsFrom("TProfile3D"):
                        continue

                    xlabel = hist.GetXaxis().GetTitle()
                    ylabel = hist.GetYaxis().GetTitle()

                    (b_min, b_max) = (hist.GetMinimumBin(), hist.GetMaximumBin())
                    (x_min, x_max) = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
                    (y_min, y_max) = (hist.GetBinContent(b_min), hist.GetBinContent(b_max))
                    y_max = 1.05*y_max

                    iPos = 33

                    # Move CMS logo/text out of frame so it does not get covered by the plot
                    if hist.InheritsFrom("TH2D") or hist.InheritsFrom("TProfile2D"):
                        iPos = 0

                    CMS.SetExtraText("Private")
                    CMS.SetEnergy("13.6")
                    CMS.SetLumi("")
                    canv = CMS.cmsCanvas('', x_min, x_max, y_min, y_max, xlabel, ylabel, square = True, extraSpace=0.06, iPos=iPos)

                    if config.has_section(hist_name):
                        logx = int(config[hist_name]["logx"])
                        logy = int(config[hist_name]["logy"])

                        if config[hist_name]["xlim"] != "":
                            (x_min, x_max) = tuple(map(float, config[hist_name]["xlim"].split(" ")))

                        if config[hist_name]["ylim"] != "":
                            (y_min, y_max) = tuple(map(float, config[hist_name]["ylim"].split(" ")))

                        if config[hist_name]["xlabel"] != "":
                            xlabel = config[hist_name]["xlabel"]

                        if config[hist_name]["ylabel"] != "":
                            xlabel = config[hist_name]["ylabel"]

                        if config[hist_name]["iPos"] != "":
                            iPos = int(config[hist_name]["iPos"])

                        canv = CMS.cmsCanvas('', x_min, x_max, y_min, y_max, xlabel, ylabel, square = True, extraSpace=0.06, iPos=0)

                        if logx == 1:
                            canv.SetLogx()

                        if logy == 1:
                            canv.SetLogy()

                        if config[hist_name]["xlabelsize"] != "":
                            xlabelsize = float(config[hist_name]["xlabelsize"])
                            CMS.GetcmsCanvasHist(canv).GetXaxis().SetLabelSize(xlabelsize)

                        if config[hist_name]["ylabelsize"] != "":
                            ylabelsize = float(config[hist_name]["ylabelsize"])
                            CMS.GetcmsCanvasHist(canv).GetYaxis().SetLabelSize(ylabelsize)

                    CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleOffset(1.5)
                    CMS.GetcmsCanvasHist(canv).GetXaxis().SetTitleOffset(0.9)

                    marker = hist.GetMarkerStyle()
                    msize = hist.GetMarkerSize()
                    mcolor = hist.GetMarkerColor()
                    lstyle = hist.GetLineStyle()
                    lwidth = hist.GetLineWidth()
                    lcolor = hist.GetLineColor()
                    fstyle = hist.GetFillStyle()
                    fcolor = hist.GetFillColor()
                    CMS.cmsDraw(hist, "", marker, msize, mcolor, lstyle, lwidth, lcolor, fstyle, fcolor)
                    canv.Draw()

                    CMS.SaveCanvas(canv, "{}/{}/{}/{}/{}.pdf".format(output_path, trigger, system_name, method_name, hist_name))

if __name__ == "__main__":
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

    config_file = args.config
    config = read_config_file(config_file)



    print("Producing plots...")
    start = time.time()
    root_files = [ROOT.TFile.Open(path) for path in files]
    for root_file in root_files:
        produce_plots(root_file, output_path, config, trigger_list=trigger_list)

    for root_file in root_files:
        root_file.Close()
    print("Finished producing plots (execution time {:.3} s)".format(time.time()-start))