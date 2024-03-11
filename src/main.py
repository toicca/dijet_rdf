import ROOT
from typing import List
import os
import socket
import argparse
from SampleAnalyzers.Dijet import DijetAnalyzer as dijet
from SampleAnalyzers.Multijet import MultijetAnalyzer as multijet
from RDFAnalyzer import JEC_corrections
from RDFHelpers import findFiles, readTriggerList

from filewriter import FileWriter
import configparser

nThreads = 4

RDataFrame = ROOT.RDataFrame

def read_config_file(config_file: str) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def parse_arguments():
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('--config', type=str, help='Path to the config file', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    config = read_config_file(args.config)
    
    print("Reading config file")
    filelist = findFiles(config["General"]["filepath"])
    nFiles = int(config["General"]["number_of_files"])
    if nFiles == -1:
        nFiles = len(filelist)
    triggerlist = readTriggerList(config["General"]["triggerpath"])
    json_file = config["General"]["golden_json"]
    nThreads = int(config["Distributed"]["nThreads"])
    ROOT.EnableImplicitMT(nThreads)
    
    print("Creating analysis object")

    corrections = JEC_corrections("", config["JEC"]["L2Relative"], "")
    dijet_analysis = dijet(filelist, triggerlist, json_file, nFiles=nFiles, JEC=corrections, nThreads=nThreads)
    multijet_analysis = multijet(filelist, triggerlist, json_file, nFiles=nFiles, JEC=corrections, nThreads=nThreads)

    dijet_analysis.do_inclusive()
    dijet_analysis.do_DB()
    dijet_analysis.do_MPF()
    dijet_analysis.run_histograms()
    hists = dijet_analysis.get_histograms()
    
    filewriter = FileWriter(config["General"]["output"])
    for trigger, histograms in hists.items():
        filewriter.write_trigger(trigger, histograms)
    filewriter.close()
    print("Done")
