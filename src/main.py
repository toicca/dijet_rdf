import ROOT
from typing import List
import os
import socket
import argparse
from SampleAnalyzers.Dijet import DijetAnalyzer as dijet
from SampleAnalyzers.Multijet import MultijetAnalyzer as multijet
from RDFAnalyzer import JEC_corrections
from RDFHelpers import findFiles, readTriggerList, parse_arguments

from filewriter import FileWriter
import configparser

nThreads = 4

RDataFrame = ROOT.RDataFrame

if __name__ == "__main__":
    # Tässä koko blockissa on jotain väärää
    args = parse_arguments()
    if args.filepath:
        filelist = findFiles(args.filepath)
    elif args.filelist != None and len(args.filelist) > 0:
        filelist = [file.strip() for file in args.filelist.split(",")]
    else:
        raise ValueError("No file list provided")
    if args.triggerpath:
        triggerlist = readTriggerList(args.triggerpath)
    elif args.triggerlist != None and len(args.triggerlist) > 0:
        triggerlist = args.triggerlist
    else:
        triggerlist = []
    nFiles = args.number_of_files
    if not nFiles:
        nFiles = len(filelist)
    elif nFiles < 0:
        nFiles = len(filelist)

    output_path = args.output_path
    run_id = args.run_id
    is_mc = args.is_MC
    is_local = args.is_local
    json_file = args.golden_json
    jetvetomap = args.jetvetomap
    L1FastJet = args.L1FastJet
    L2Relative = args.L2Relative
    L2L3Residual = args.L2L3Residual
    JER = args.JER
    JER_SF = args.JER_SF
    nThreads = args.nThreads
    verbosity = args.verbosity
    progress_bar = args.progress_bar
    cutflow_report = args.cutflow_report
    
    ROOT.EnableImplicitMT(nThreads)
    
    print("Creating analysis object")
    print(L2Relative, args.L2Relative, args.L2Relative)

    corrections = JEC_corrections(L1FastJet, L2Relative, L2L3Residual, JER, JER_SF)
    dijet_analysis = dijet(filelist, triggerlist, json_file, nFiles=nFiles, JEC=corrections, nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local)

    dijet_analysis.do_inclusive()
    dijet_analysis.do_PFComposition()
    dijet_analysis.do_DB()
    dijet_analysis.do_MPF()
    dijet_analysis.do_RunsAndLumis()
    if is_mc:
        dijet_analysis.do_MC()
    dijet_analysis.run_histograms()
    hists = dijet_analysis.get_histograms()
    
    filewriter = FileWriter(output_path + "/dijet_" + run_id + ".root")
    for trigger, histograms in hists.items():
        filewriter.write_trigger(trigger, histograms)
    filewriter.close()
    print("Done")
