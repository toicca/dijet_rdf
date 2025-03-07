import argparse
import ROOT

import find_json
import find_newest
import find_range
import histograms
import produce_ratio
import produce_responses
import produce_time_evolution
import produce_vetomaps
from plotting import produce_plots
import skim

class ProcessingState:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='JEC4PROMPT Toolkit: \
            https://github.com/toicca/dijet_rdf/tree/main')
        self.subparsers = self.parser.add_subparsers(dest='subparser_name')
        self.args = None
        self.valfuncs = {}
        self.commands = {}

    def init_state(self):
        skim.update_state(self)
        histograms.update_state(self)
        # histograms.add_hist_parser(self.subparsers)
        # find_json.add_find_json_parser(self.subparsers)
        # find_newest.add_find_newest_parser(self.subparsers)
        # find_range.add_find_range_parser(self.subparsers)
        # produce_plots.add_plots_parser(self.subparsers)
        # produce_vetomaps.add_vetomaps_parser(self.subparsers)
        
        self.args = self.parser.parse_args()
        self.valfuncs[self.args.subparser_name](self.args)

    def execute_command(self):
        self.commands[self.args.subparser_name](self)

if __name__ == "__main__":
    state = ProcessingState()
    state.init_state()
    state.execute_command()
