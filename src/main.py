import ROOT
from SampleAnalyzers.Dijet import DijetAnalyzer as dijet
from SampleAnalyzers.Multijet import MultijetAnalyzer as multijet
from SampleAnalyzers.MultijetM import MultijetMAnalyzer as multijetM
from SampleAnalyzers.DijetM import DijetMAnalyzer as dijetM
from RDFAnalyzer import RDFAnalyzer
from RDFHelpers import parse_arguments

 
verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kDebug+10)
from filewriter import FileWriter

nThreads = 4

RDataFrame = ROOT.RDataFrame

if __name__ == "__main__":
    # Tässä koko blockissa on jotain väärää
    args = parse_arguments()
    if args.filepath:
        filelist = args.filepath
    elif args.filelist:
        filelist = args.filelist
    else:
        raise ValueError("No file list provided")

    if args.triggerpath:
        triggerlist = args.triggerpath
        triggers = {}
        for trigger in triggerlist:
            triggers[trigger] = trigger
    elif args.triggerlist:
        triggerlist = args.triggerlist
        triggers = {}
        for trigger in triggerlist:
            triggers[trigger] = trigger
    elif args.triggerconfig:
        triggers = args.triggerconfig

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

    ROOT.EnableImplicitMT(nThreads)

    print("Creating analysis object " + ("(MC)" if is_mc else "(Data)"))
    
    standard_analysis = RDFAnalyzer(filelist, triggers, json_file, nFiles=nFiles, JEC=correction_dict, \
                                    nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local, run_raw=run_raw, \
                                    selection_only=selection_only, header_dir=header_dir)
    
    dijet_analysis = dijetM(filelist, triggers, json_file, nFiles=nFiles, JEC=correction_dict, \
                           nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local, run_raw=run_raw, \
                           selection_only=selection_only, header_dir=header_dir)
    
    # multijet_analysis = multijet(filelist, triggers, json_file, nFiles=nFiles, JEC=correction_dict, \
                                # nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local, run_raw=run_raw, \
                                # selection_only=selection_only, header_dir=header_dir)

    multijet_analysis = multijetM(filelist, triggers, json_file, nFiles=nFiles, JEC=correction_dict, \
                                nThreads=nThreads, progress_bar=progress_bar, isMC=is_mc, local=is_local, run_raw=run_raw, \
                                selection_only=selection_only, header_dir=header_dir)
    
    # standard_analysis.do_inclusive()
    # standard_analysis.do_inclusive_control()
    # standard_analysis.do_PFComposition()
    # standard_analysis.do_RunsAndLumis()
    
    if is_mc:
        standard_analysis.do_MC()
    
    dijet_analysis.do_DB()
    dijet_analysis.do_sample_control()
    dijet_analysis.do_MPF()
    multijet_analysis.do_sample_control()
    multijet_analysis.do_DB()
    multijet_analysis.do_MPF()
    
    standard_analysis.run_histograms()
    dijet_analysis.run_histograms()
    multijet_analysis.run_histograms()
    
    hists1 = standard_analysis.get_histograms()
    hists2 = dijet_analysis.get_histograms()
    hists3 = multijet_analysis.get_histograms()
    
    filewriter = FileWriter(output_path + "/multisample_" + run_id + ".root", standard_analysis.trigger_list, cut_hist_names)
    filewriter.write_samples([standard_analysis, dijet_analysis, multijet_analysis])
    filewriter.close()
    print("Done")
