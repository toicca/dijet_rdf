import ROOT
from SampleAnalyzers.Dijet import DijetAnalyzer as dijet
from SampleAnalyzers.Multijet import MultijetAnalyzer as multijet
from RDFAnalyzer import RDFAnalyzer
from RDFHelpers import parse_arguments
 
import pandas as pd
import numpy as np

verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kDebug+10)
from filewriter import FileWriter

import analysis
import find_json
import find_newest
import find_range
import produce_ratio
import produce_responses
import produce_time_evolution
import produce_vetomaps
from plotting import produce_plots

nThreads = 4

RDataFrame = ROOT.RDataFrame

if __name__ == "__main__":
    # Tässä koko blockissa on jotain väärää
    args = parse_arguments()

    command = args.subparser_name
    
    # One day when Python version >= 3.10.0
    # is used implement match here instead.
    if command == "analysis":
        analysis.run(args)
    elif command == "find_json":
        find_json.run(args)
    elif command == "find_newest":
        find_newest.run(args)
    elif command == "find_range":
        find_range.run(args)
    elif command == "produce_ratio":
        produce_ratio.run(args)
    elif command ==  "produce_responses":
        produce_responses.run(args)
    elif command == "produce_time_evolution":
        produce_time_evolution.run(args)
    elif command == "produce_vetomaps":
        produce_vetomaps.run(args)
    elif command == "produce_plots":
        produce_plots.run(args)
