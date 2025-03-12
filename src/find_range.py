import ROOT
import argparse

from utils.processing_utils import file_read_lines

def update_state(state):
    add_find_range_parser(state.subparsers)
    state.valfuncs['find_range'] = validate_args
    state.commands['find_range'] = run

def add_find_range_parser(subparsers):
    find_range_parser = subparsers.add_parser('find_range', help="Find run range of \
                                                given input files")
    find_range_files = find_range_parser.add_mutually_exclusive_group(required=True)
    find_range_files.add_argument("--filelist", type=str, help="Comma separated list of \
            input files")
    find_range_files.add_argument('-fp', '--filepaths', type=str, help='Comma separated list of \
            text files containing input files (one input file per line).')
    find_range_parser.add_argument("-loc", "--is_local", action="store_true", help='Run locally. \
            If not set will append root://cms-xrd-global.cern.ch/ \
            to the start of file names')
    find_range_parser.add_argument("--for_brilcalc", action="store_true", help='Prints the range \
            in a form compatible with the brilcalc command line tool')
    find_range_parser.add_argument("--nThreads", type=int, help="Number of threads to be used \
            for multithreading")
    find_range_parser.add_argument("--progress_bar", action="store_true", help="Show progress bar")

def validate_args(args):
    pass

def find_run_range(rdf):
    return int(rdf.Min("run").GetValue()), int(rdf.Max("run").GetValue())

def run(args):
    # Shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    # Split the file list and trigger list if they are given as a string
    if args.filelist:
        files= args.filelist.split(",")
    elif args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        for path in paths:
            files.extend(file_read_lines(path, find_ROOT=True))
    else:
        raise ValueError("No file list provided")

    chain = ROOT.TChain("Events")
    for file in files:
        if args.is_local:
            chain.Add(file)
        else:
            chain.Add(f"root://cms-xrd-global.cern.ch/{file}")

    rdf = ROOT.RDataFrame(chain)
    if args.progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(rdf)

    min_run, max_run = find_run_range(rdf)
    if args.for_brilcalc:
        print(f"--begin {min_run} --end {max_run}");
    else:
        print(f"{min_run},{max_run}")
