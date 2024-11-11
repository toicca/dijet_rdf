import ROOT
import argparse

from RDFHelpers import file_read_lines

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
    elif args.filepath:
        files = file_read_lines(args.filepath, find_ROOT=True)
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
