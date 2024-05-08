import ROOT
import ctypes

from RDFHelpers import file_read_lines

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run range finder for dijet_rdf: https://github.com/toicca/dijet_rdf")

    files = parser.add_mutually_exclusive_group(required=True)
    files.add_argument("--filelist", type=str, help="Comma separated list of input files")
    files.add_argument("--filepath", type=str, help="Path to a root file containing a list of input files")

    parser.add_argument("-loc", "--is_local", action="store_true", help='Run locally. If not set will append root://cms-xrd-global.cern.ch/ to the start of file names')

    args = parser.parse_args()
    
    return args

def find_run_range(rdf):
    return int(rdf.Min("run").GetValue()), int(rdf.Max("run").GetValue())

if __name__ == "__main__":
    args = parse_arguments()
    if args.filepath:
        files = file_read_lines(args.filepath)
    else:
        files = [s.strip() for s in args.filelist.split(',')]

    chain = ROOT.TChain("Events")
    for file in files:
        if args.is_local:
            chain.Add(f"root://cms-xrd-global.cern.ch/{file}")
        else:
            chain.Add(file)

    rdf = ROOT.RDataFrame(chain)
    min_run, max_run = find_run_range(rdf)
    print(min_run, max_run)
