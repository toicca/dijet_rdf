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
    files.add_argument("--filepath", type=str, help="Path to a text file containing a list of output files produced by dijet_rdf")

    parser.add_argument("--out", required=True, type=str, help="Output path")

    parser.add_argument("--config", type=str, default="", help="Path to config file")

    parser.add_argument("--all", action="store_true", help="Produce all plots in given .root files")

    args = parser.parse_args()
    
    return args

def read_config_file(file: str) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.read(file)

    return config

def produce_plots_all(file, output_path):
    subfolder_name = file.GetName().split("/")[-1].replace(".root", "").replace(".", "_")
    for trigger_key in file.GetListOfKeys():
        trigger_name = trigger_key.GetName()
        systems = file.Get(trigger_name)
        for system_key in systems.GetListOfKeys():
            system_name = system_key.GetName()
            methods = systems.Get(system_name)
            for method_key in methods.GetListOfKeys():
                method_name = method_key.GetName()

                pathlib.Path("{}/{}/{}/{}/{}".format(output_path, subfolder_name, trigger_name, system_name, method_name)).mkdir(exist_ok=True, parents=True)

                hists = methods.Get(method_name)
                for hist_key in hists.GetListOfKeys():
                    hist_name = hist_key.GetName()
                    hist = hists.Get(hist_name)

                    # ignore 3D histograms and profiles for now
                    if hist.InheritsFrom("TH3D") or hist.InheritsFrom("TProfile3D"):
                        continue

                    xtitle = hist.GetXaxis().GetTitle()
                    ytitle = hist.GetYaxis().GetTitle()

                    (b_min, b_max) = (hist.GetMinimumBin(), hist.GetMaximumBin())
                    (x_min, x_max) = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
                    (y_min, y_max) = (hist.GetBinContent(b_min), hist.GetBinContent(b_max))
                    y_max = 1.05*y_max
                    iPos = 33

                    # Move CMS logo/text out of frame so it does not get covered by the plot
                    if hist.InheritsFrom("TH2D") or hist.InheritsFrom("TProfile2D"):
                        (x_min, x_max) = (hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
                        (y_min, y_max) = (hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax())
                        iPos = 0

                    CMS.SetExtraText("Preliminary")
                    CMS.SetEnergy("13.6")
                    CMS.SetLumi("")
                    canv = CMS.cmsCanvas('', x_min, x_max, y_min, y_max, xtitle, ytitle, square = True, extraSpace=0.0, iPos=iPos)

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

                    CMS.SaveCanvas(canv, "{}/{}/{}/{}/{}/{}.pdf".format(output_path, subfolder_name, trigger_name, system_name, method_name, hist_name))

def produce_plots_from_config(file, output_path, config, plots):
    subfolder_name = file.GetName().split("/")[-1].replace(".root", "").replace(".", "_")
    for plot_name in plots:
        if config[plot_name]["triggers"] == "":
            trigger_keys = file.GetListOfKeys()
            triggers = [tkey.GetName() for tkey in trigger_keys]
        else:
            triggers = [s.strip() for s in config[plot_name]["triggers"].split(",")]

        hists_config = [s.strip() for s in config[plot_name]["hists"].split(",")]

        for trigger in triggers:
            systems = file.Get(trigger)
            collected_plots = []
            for system_key in systems.GetListOfKeys():
                system_name = system_key.GetName()
                methods = systems.Get(system_name)
                for method_key in methods.GetListOfKeys():
                    method_name = method_key.GetName()
                    hists = methods.Get(method_name)
                    for hist_key in hists.GetListOfKeys():
                        hist_name = hist_key.GetName()
                        method_hist_str = f"{method_name}/{hist_name}"
                        system_method_hist_str = f"{system_name}/{method_name}/{hist_name}"
                        #print(system_name, method_name, hist_name)
                        if hist_name in hists_config or method_hist_str in hists_config or system_method_hist_str in hists_config:
                            if config[plot_name]["projection"] == "x":
                                collected_plots.append(hists.Get(hist_name).ProjectionX())
                            elif config[plot_name]["projection"] == "y":
                                collected_plots.append(hists.Get(hist_name).ProjectionY())
                            elif config[plot_name]["projection"] == "z":
                                collected_plots.append(hists.Get(hist_name).ProjectionZ())
                            else:
                                collected_plots.append(hists.Get(hist_name))
            
            if len(collected_plots) == 0: continue

            all_1D = all(p.InheritsFrom("TH1D") or p.InheritsFrom("TProfile1D") for p in collected_plots)
            all_2D = all(p.InheritsFrom("TH2D") or p.InheritsFrom("TProfile2D") for p in collected_plots)

            if not (all_1D or all_2D):
                print(f"Skipping plotting of [{plot_name}], histogram list given in the config file contains both 1D and 2D histograms")

            pathlib.Path("{}/{}/{}".format(output_path, subfolder_name, trigger)).mkdir(exist_ok=True, parents=True)

            xtitle = collected_plots[0].GetXaxis().GetTitle()
            ytitle = collected_plots[0].GetYaxis().GetTitle()

            (x_min, x_max) = (min([p.GetXaxis().GetXmin() for p in collected_plots]), max([p.GetXaxis().GetXmax() for p in collected_plots]))
            (y_min, y_max) = (min([p.GetYaxis().GetXmin() for p in collected_plots]), max([p.GetYaxis().GetXmax() for p in collected_plots]))

            if all_1D:
                (x_min, x_max) = (min([p.GetXaxis().GetXmin() for p in collected_plots]), max([p.GetXaxis().GetXmax() for p in collected_plots]))
                (y_min, y_max) = (min([p.GetBinContent(p.GetMinimumBin()) for p in collected_plots]), max([p.GetBinContent(p.GetMaximumBin()) for p in collected_plots]))
                y_max = 1.05*y_max

            iPos = 33
            if all_2D:
                iPos = 0

            extraSpace = 0.0

            CMS.SetExtraText("Preliminary")
            CMS.SetEnergy("13.6")
            if config[plot_name]["energy"] != "":
                CMS.SetEnergy(config[plot_name]["energy"])
            
            CMS.SetLumi("")
            logx = int(config[plot_name]["logx"])
            logy = int(config[plot_name]["logy"])

            if config[plot_name]["xlim"] != "":
                (x_min, x_max) = tuple(map(float, config[plot_name]["xlim"].split(" ")))

            if config[plot_name]["ylim"] != "":
                (y_min, y_max) = tuple(map(float, config[plot_name]["ylim"].split(" ")))

            if config[plot_name]["xtitle"] != "":
                xtitle = config[plot_name]["xtitle"]

            if config[plot_name]["ytitle"] != "":
                ytitle = config[plot_name]["ytitle"]

            if config[plot_name]["iPos"] != "":
                iPos = int(config[plot_name]["iPos"])

            if config[plot_name]["extraSpace"] != "":
                extraSpace = float(config[plot_name]["extraSpace"])

            canv = CMS.cmsCanvas('', x_min, x_max, y_min, y_max, xtitle, ytitle, square=True, extraSpace=extraSpace, iPos=iPos)
            if logx == 1:
                canv.SetLogx()
            if logy == 1:
                canv.SetLogy()

            if config[plot_name]["xlabelsize"] != "":
                xlabelsize = float(config[plot_name]["xlabelsize"])
                CMS.GetcmsCanvasHist(canv).GetXaxis().SetLabelSize(xlabelsize)

            if config[plot_name]["ylabelsize"] != "":
                ylabelsize = float(config[plot_name]["ylabelsize"])
                CMS.GetcmsCanvasHist(canv).GetYaxis().SetLabelSize(ylabelsize)

            if config[plot_name]["xtitlesize"] != "":
                xtitlesize = float(config[plot_name]["xtitlesize"])
                CMS.GetcmsCanvasHist(canv).GetXaxis().SetTitleSize(xtitlesize)

            if config[plot_name]["ytitlesize"] != "":
                ytitlesize = float(config[plot_name]["ytitlesize"])
                CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleSize(ytitlesize)

            if config[plot_name]["xtitleoffset"] != "":
                xtitleoffset = float(config[plot_name]["xtitleoffset"])
                CMS.GetcmsCanvasHist(canv).GetXaxis().SetTitleOffset(xtitleoffset)

            if config[plot_name]["ytitleoffset"] != "":
                ytitleoffset = float(config[plot_name]["ytitleoffset"])
                CMS.GetcmsCanvasHist(canv).GetYaxis().SetTitleOffset(ytitleoffset)

            markers = []
            if config[plot_name]["markers"] != "":
                markers = [s.strip() for s in config[plot_name]["markers"].split(",")]
            
            colors = []
            if config[plot_name]["markers"] != "":
                colors = [s.strip() for s in config[plot_name]["colors"].split(",")]

            leg = None
            entries = []
            if config[plot_name]["legend"] != "" and config[plot_name]["legendPos"] != "":
                entries = [s.strip() for s in config[plot_name]["legend"].split("@")]
                pos = [float(s.strip()) for s in config[plot_name]["legendPos"].split(",")]
                leg = CMS.cmsLeg(pos[0], pos[1], pos[2], pos[3])

            for (i, plot) in enumerate(collected_plots):

                if leg is not None and len(entries) > 0:
                    leg.AddEntry(plot, entries[i])

                if (len(markers) > 0):
                    marker = int(markers[i])
                    fstyle = int(markers[i])
                else:
                    marker = plot.GetMarkerStyle()
                    fstyle = plot.GetFillStyle()

                if (len(colors) > 0):
                    mcolor = int(colors[i])
                    lcolor = int(colors[i])
                    fcolor = int(colors[i])
                else:
                    mcolor = plot.GetMarkerColor()
                    lcolor = plot.GetLineColor()
                    fcolor = plot.GetFillColor()
                
                lstyle = plot.GetLineStyle()
                msize = plot.GetMarkerSize()
                lwidth = plot.GetLineWidth()
                CMS.cmsDraw(plot, "", marker=marker, msize=msize, mcolor=mcolor, lstyle=lstyle, lwidth=lwidth, lcolor=lcolor, fstyle=fstyle, fcolor=fcolor)

            CMS.SaveCanvas(canv, "{}/{}/{}/{}.pdf".format(output_path, subfolder_name, trigger, plot_name))

if __name__ == "__main__":
    args = parse_arguments()
    
    files: List[str] = []

    if args.filepath:
        files = file_read_lines(args.filepath)
    else:
        files = [s.strip() for s in args.filelist.split(',')]

    output_path = args.out

    config_file = args.config
    config = read_config_file(config_file)

    plots = config.sections()
    print("Producing plots...")
    start = time.time()
    root_files = [ROOT.TFile.Open(path) for path in files]
    for root_file in root_files:
        produce_plots_from_config(root_file, output_path, config, plots)
        if args.all:
            produce_plots_all(root_file, output_path)

    for root_file in root_files:
        root_file.Close()
    print("Finished producing plots (execution time {:.3} s)".format(time.time()-start))
