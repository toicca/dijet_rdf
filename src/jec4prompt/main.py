import argparse
import logging
from pathlib import Path

import jec4prompt.find_json as find_json
import jec4prompt.find_newest as find_newest
import jec4prompt.find_range as find_range
import jec4prompt.histograms as histograms
import jec4prompt.produce_ratio as produce_ratio
import jec4prompt.produce_responses as produce_responses
import jec4prompt.produce_time_evolution as produce_time_evolution
import jec4prompt.produce_vetomaps as produce_vetomaps
import jec4prompt.skim as skim
from jec4prompt.plotting import produce_plots


class ProcessingState:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description="JEC4PROMPT Toolkit: \
            https://github.com/toicca/dijet_rdf/tree/main"
        )
        self.subparsers = self.parser.add_subparsers(dest="subparser_name")
        self.parser.add_argument(
            "--log",
            type=str,
            default="INFO",
            help="Logging level",
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        )
        self.parser.add_argument(
            "--tag", type=str, help="Processing tag for the output file", default=""
        )
        self.parser.add_argument(
            "--redirector",
            type=str,
            help="Redirector for the input files",
            default="root://xrootd-cms.infn.it//",
        )
        self.args = None
        self.channels = [
            "photonjet",
            "dijet",
            "multijet",
            "zmm",
            "zee",
            "zjet",
            "empty",
        ]
        self.valfuncs = {}
        self.commands = {}
        self.module_dir = Path(__file__).resolve().parent
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

    def init_state(self):
        self.logger.info("Initializing state")
        skim.update_state(self)
        histograms.update_state(self)
        find_json.update_state(self)
        find_newest.update_state(self)
        find_range.update_state(self)
        produce_ratio.update_state(self)
        produce_responses.update_state(self)
        produce_time_evolution.update_state(self)
        produce_plots.update_state(self)
        produce_vetomaps.update_state(self)

        self.args = self.parser.parse_args()
        self.logger.info(f"Parsed arguments: {self.args}")
        self.logger.info(f"Subparser name: {self.args.subparser_name}")
        self.logger.info(f"Logging level: {self.args.log}")
        self.logger.setLevel(self.args.log)
        
        # Check that there's a subparser
        if self.args.subparser_name is None:
            self.parser.print_help()
            self.parser.exit()

        self.valfuncs[self.args.subparser_name](self.args)

    def execute_command(self):
        self.logger.info(f"Executing command: {self.args.subparser_name}")
        self.commands[self.args.subparser_name](self)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


def main():
    state = ProcessingState()
    state.init_state()
    original_state = state
    state.execute_command()
    if state != original_state:
        original_state.logger.error("State has changed after execution")


if __name__ == "__main__":
    main()
