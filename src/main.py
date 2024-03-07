import ROOT
from typing import List
import os
from distributed import Client
from dask_lxplus import CernCluster
import socket
import argparse
from dijet import dijet, JEC_corrections

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

def findFiles(file : str) -> List[str]:
    with open(file) as f:
        return [line.strip() for line in f.readlines() if line.strip().endswith(".root")]
    
def readTriggerList(file : str) -> List[str]:
    with open(file) as f:
        return [line.strip() for line in f.readlines()]

def create_connection(cores):
    n_port = 8786
    cluster = CernCluster(
            cores=cores,
            memory='2000MB',
            disk='1000MB',
            death_timeout = '60',
            lcg = True,
            nanny = False,
            container_runtime = "none",
            log_directory = "/eos/user/d/dnasman/condor/log",
            scheduler_options={
                'port': n_port,
                'host': socket.gethostname(),
                },
            job_extra={
                '+JobFlavour': '"tomorrow"',
                },
            extra = ['--worker-port 10000:10100']
            )
    try:
        client = Client(cluster, timeout='2s')
    except TimeoutError:
        pass
    return client, cluster

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

    dijet_analysis.do_inclusive()
    dijet_analysis.do_DB()
    dijet_analysis.do_DB(system="multijet")
    dijet_analysis.do_MPF()
    dijet_analysis.do_MPF(system="multijet")
    dijet_analysis.run_histograms()
    hists = dijet_analysis.get_histograms()
    
    filewriter = FileWriter(config["General"]["output"])
    for trigger, histograms in hists.items():
        filewriter.write_trigger(trigger, histograms)
    filewriter.close()
    # Rest of your code goes here
    print("Hello World!")