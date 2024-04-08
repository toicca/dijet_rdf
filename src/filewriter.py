import ROOT
import os
from typing import List, Dict

class FileWriter:
    def __init__(self, output_file : str, triggers : List[str] = [], samples : Dict[str, List] = {}):
        self.triggers = triggers
        self.samples = samples
        
        self.output_file = output_file
        if not os.path.exists(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        self.output = ROOT.TFile(output_file, "RECREATE")
        self.output.cd()
        
    def write_trigger(self, trigger : str, histograms : list):
        self.output.mkdir(trigger)
        self.output.cd(trigger)
        for hist in histograms:
            hist_name = hist.GetName()
            file_name = hist_name.split("_")[0]
            if not self.output.GetDirectory(trigger + "/" + file_name):
                self.output.mkdir(trigger + "/" + file_name)
            self.output.cd(trigger + "/" + file_name)
            hist.Write()
            self.output.cd("..")
            
        self.output.cd()
        
    def write_full(self, sample : str, trigger : str, histograms : list):
        self.write_trigger(sample + "/" + trigger, histograms)
        
    def close(self):
        self.output.Close()
        