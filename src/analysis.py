import ROOT
from SampleAnalyzers.Dijet import DijetAnalyzer as dijet
from SampleAnalyzers.Multijet import MultijetAnalyzer as multijet
from RDFAnalyzer import RDFAnalyzer
from RDFHelpers import parse_arguments, read_correction_config, read_trigger_config, file_read_lines
 
import pandas as pd
import numpy as np

verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kDebug+10)
from filewriter import FileWriter

nThreads = 4

RDataFrame = ROOT.RDataFrame

def run(args):
    # If config file is set, override all arguments with ones set in there
    if args.config:
        config = read_config_file(args.config)
        for arg, value in config["GENERAL"].items():
            # Do a type conversion for the option
            if arg == "number_of_files" or arg == "nThreads" \
            or arg == "verbosity" or arg == "is_local" \
            or arg == "is_mc" or arg == "progress_bar" \
            or arg == "cutflow_report"  or arg == "cut_histogram_names" \
            or arg == "selection_only" or arg == "trigger_details":
                if value == "":
                    setattr(args, arg, 0)
                else:
                    setattr(args, arg, int(value))
            else:
                setattr(args, arg, value)
                
    # Split the file list and trigger list if they are given as a string
    if args.filelist:
        filelist= args.filelist.split(",")
    elif args.filepath:
        filelist = file_read_lines(args.filepath, find_ROOT=True)
    else:
        raise ValueError("No file list provided")
        
    if args.triggerlist:
        triggerlist = args.triggerlist.split(",")
    elif args.triggerpath:
        triggerlist = file_read_lines(args.triggerpath)
    elif args.trigger_config:
        triggerlist = read_trigger_config(args.trigger_config)

    triggers = {}
    for trigger in triggerlist:
        triggers[trigger] = trigger

    if args.correction_config:
        args.correction_config = read_correction_config(args.correction_config)

    nFiles = args.number_of_files
    if not nFiles:
        nFiles = len(filelist)
    elif nFiles < 0:
        nFiles = len(filelist)

    output_path = args.output_path
    run_id = args.run_id
    is_mc = args.is_MC
    is_local = args.is_local
    json_file = args.golden_json
    jetvetomap = args.jetvetomap
    correction_dict = args.correction_config
    nThreads = args.nThreads
    verbosity = args.verbosity
    progress_bar = args.progress_bar
    cutflow_report = args.cutflow_report
    cut_hist_names = args.cut_histogram_names
    run_raw = args.run_raw
    selection_only = args.selection_only
    header_dir = args.header_dir
    trigger_details = args.trigger_details

    if args.lumi:
        df = pd.read_csv(args.lumi, comment='#', names=["run:fill", "time", "nls", "ncms", "delivered(/fb)", "recorded(/fb)"])
        recorded_luminosity = df["recorded(/fb)"].to_numpy()
        print(f"Running on {np.sum(recorded_luminosity)} 1/fb luminosity")

    if nThreads > 0:
        ROOT.EnableImplicitMT(nThreads)
    else:
        ROOT.DisableImplicitMT()

    print("Creating analysis object " + ("(MC)" if is_mc else "(Data)"))
    
    standard_analysis = RDFAnalyzer(filelist, triggers, trigger_details=trigger_details, json_file= json_file, nFiles=nFiles, JEC=correction_dict, \
                                    nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local, run_raw=run_raw, \
                                    selection_only=selection_only, header_dir=header_dir)
    
    dijet_analysis = dijet(filelist, triggers, trigger_details=trigger_details, json_file= json_file, nFiles=nFiles, JEC=correction_dict, \
                           nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local, run_raw=run_raw, \
                           selection_only=selection_only, header_dir=header_dir)
    
    multijet_analysis = multijet(filelist, triggers, trigger_details=trigger_details, json_file= json_file, nFiles=nFiles, JEC=correction_dict, \
                                nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local, run_raw=run_raw, \
                                selection_only=selection_only, header_dir=header_dir)

    standard_analysis.do_inclusive()
    standard_analysis.do_inclusive_control()
    standard_analysis.do_PFComposition()
    standard_analysis.do_RunsAndLumis()
    
    if is_mc:
        standard_analysis.do_MC()
    
    dijet_analysis.do_DB()
    dijet_analysis.do_sample_control()
    dijet_analysis.do_MPF()
    dijet_analysis.do_Activity()
    multijet_analysis.do_sample_control()
    multijet_analysis.do_DB()
    multijet_analysis.do_MPF()
    
    standard_analysis.run_histograms()
    dijet_analysis.run_histograms()
    multijet_analysis.run_histograms()

    hists1 = standard_analysis.get_histograms()
    hists2 = dijet_analysis.get_histograms()
    hists3 = multijet_analysis.get_histograms()
    
    filewriter = FileWriter(output_path + f"/{run_id}_era{standard_analysis.era}_runs{standard_analysis.run_range[0]}to{standard_analysis.run_range[1]}.root", standard_analysis.histograms.keys(), cut_hist_names)
    filewriter.write_samples([standard_analysis, dijet_analysis, multijet_analysis])
    filewriter.close()
    print("Done")
