import os
from typing import List

import ROOT
from RDFAnalyzer import RDFAnalyzer


class FileWriter:
    def __init__(
        self, output_file: str, triggers: List[str] = None, cut_hist_names: bool = False
    ):
        self.triggers = triggers
        self.cut_hist_names = cut_hist_names

        self.output_file = output_file
        if not os.path.exists(os.path.dirname(output_file)):
            os.makedirs(os.path.dirname(output_file))
        self.output = ROOT.TFile(output_file, "RECREATE")
        self.output.cd()

    def write_trigger(self, trigger: str, histograms: list):
        self.output.mkdir(trigger)
        self.output.cd(trigger)
        for hist in histograms:
            hist_name = hist.GetName()
            file_name = hist_name.split("_")[0]
            if not self.output.GetDirectory(trigger + "/" + file_name):
                self.output.mkdir(trigger + "/" + file_name)
            self.output.cd(trigger + "/" + file_name)
            if self.cut_hist_names:
                hist_name = hist_name.split("_")[-1]
                hist.SetName(hist_name)
            hist.Write()
            self.output.cd("..")

        self.output.cd("..")

    def write_samples(self, samples: List[RDFAnalyzer], triggers: List[str] = None):

        if not triggers:
            triggers = self.triggers

        for trigger in triggers:
            for sample in samples:
                histograms = sample.get_histograms()

                for hist in histograms[trigger]:
                    hist_name = hist.GetName()
                    file_name = hist_name.split("_")[0]
                    if not self.output.GetDirectory(
                        trigger + "/" + sample.system + "/" + file_name
                    ):
                        self.output.mkdir(
                            trigger + "/" + sample.system + "/" + file_name
                        )
                    self.output.cd(trigger + "/" + sample.system + "/" + file_name)
                    if self.cut_hist_names:
                        if hist_name.split("_")[-1] == "selected":
                            hist_name = hist_name.split("_")[-2] + "_selected"
                        elif hist_name.split("_")[-1] == "all":
                            hist_name = hist_name.split("_")[-2] + "_all"
                        else:
                            hist_name = hist_name.split("_")[-1]

                    hist.SetName(hist_name)
                    hist.Write()
                    self.output.cd()

    def close(self):
        self.output.Close()
