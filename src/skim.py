import ROOT
import os
import subprocess
import argparse
import pathlib
import ctypes
import numpy as np
import pandas as pd
import time
from typing import List

from processing_utils import file_read_lines, find_site
from skimming_utils import filter_json, correct_jets, find_vetojets, get_Flags, sort_jets

weight_info = {
    "xsec" : {
        # dijet and multijet
        # QCD xsecs the same for Summer24 and Winter24
        "QCD-4Jets_HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 892.4,
        "QCD-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 25240000.0,
        "QCD-4Jets_HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8": 385.4,
        "QCD-4Jets_HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 126.5,
        "QCD-4Jets_HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 26.53,
        "QCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1958000.0,
        "QCD-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 96730.0,
        "QCD-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": 312200000.0,
        "QCD-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8": 13590.0,
        "QCD-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": 58840000.0,
        "QCD-4Jets_HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 3046.0,
        # zjet
        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": 6695.0, # Winter24
        "DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1.0, # Summer24 -> single sample so xsec doesn't contribute
        # egamma Summer24
        "GJ-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": 15240.0,
        "GJ-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": 8111.0,
        "GJ-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 7327.0,
        "GJ-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1541.0,
        "GJ-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 167.6,
        "GJ-4Jets_HT-600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 54.39,
        # egamma Winter24
        "GJ-4Jets_Bin-HT-40to100-PTG-10to100_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 123200.0,
        "GJ-4Jets_Bin-HT-100to200-PTG-10to100_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 32190.0,
        "GJ-4Jets_Bin-HT-200to400-PTG-10to100_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 5514.0,
        "GJ-4Jets_Bin-HT-400to600-PTG-10to100_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 483.8,
        "GJ-4Jets_Bin-HT-600to1000-PTG-10to100_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 117.4,
        "GJ-4Jets_Bin-HT-1000-PTG-10to100_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 15.11,
        "GJ-4Jets_Bin-HT-40to200-PTG-100to200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 557.0,
        "GJ-4Jets_Bin-HT-200to400-PTG-100to200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 202.4,
        "GJ-4Jets_Bin-HT-400to600-PTG-100to200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 29.95,
        "GJ-4Jets_Bin-HT-600to1000-PTG-100to200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 9.646,
        "GJ-4Jets_Bin-HT-1000-PTG-100to200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1.632,
        "GJ-4Jets_Bin-HT-40to400-PTG-200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 43.92,
        "GJ-4Jets_Bin-HT-400to600-PTG-200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 11.77,
        "GJ-4Jets_Bin-HT-600to1000-PTG-200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 4.743,
        "GJ-4Jets_Bin-HT-1000-PTG-200_Par-dRGJ-0p25_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1.018
    },
# TODO
#    "nGenEvents" : {
#        # dijet and multijet
#        "QCD-4Jets_HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 2895970,
#        "QCD-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 5629540,
#        "QCD-4Jets_HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19537600,
#        "QCD-4Jets_HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 17527100,
#        "QCD-4Jets_HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 9212540,
#        "QCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": 18647200,
#        "QCD-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19101200,
#        "QCD-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19282700,
#        "QCD-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8": 19122400,
#        "QCD-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1054540,
#        "QCD-4Jets_HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8": 18625600,
#        # zjet
#        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": 74301700
#        # egamma
#        "GJ-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???,
#        "GJ-4Jets_HT-600_TuneCP5_13p6TeV_madgraphMLM-pythia8": ???
#    }
}

jet_columns = [
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId",
    "Jet_area", "Jet_nConstituents", "Jet_nElectrons", "Jet_nMuons",
    "Jet_chEmEF", "Jet_chHEF",
    "Jet_neEmEF", "Jet_neHEF",
    "Jet_hfEmEF", "Jet_hfHEF",
    "Jet_muEF",
    "Jet_neMultiplicity", "Jet_chMultiplicity",
    "Jet_rawFactor",
    "Jet_btagPNetQvG"
]


def init_TnP(rdf, channel):
    """
    Initialize the tag and probe variables for the analysis
    """
    if channel == "dijet":
        from selections.dijet import init_dijet as init_selection
    elif channel == "zjet":
        from selections.zmm import init_zmm as init_selection
    elif channel == "egamma":
        from selections.photonjet import init_photonjet as init_selection
    elif channel == "multijet":
        from selections.multijet import init_multijet as init_selection
    else:
        from selections.empty import init_empty as init_selection

    rdf = init_selection(rdf, jet_columns)

    rdf = rdf.Define("Probe_rawPt", "(1.0 - Probe_rawFactor) * Probe_pt")
    rdf = rdf.Filter("Activity_idx_temp >= 0 ? Jet_pt[Activity_idx_temp] / ((Probe_pt + Tag_pt)*0.5) < 1.0 : 1", "Activity jet pT fraction < 1.0")

    # Label non-flat branches as _temp to drop them later
    rdf = (rdf.Define("Tag_fourVec_temp", "ROOT::Math::PtEtaPhiMVector(Tag_pt, Tag_eta, Tag_phi, Tag_mass)")
            .Define("Probe_fourVec_temp", "ROOT::Math::PtEtaPhiMVector(Probe_pt, Probe_eta, \
                    Probe_phi, Probe_mass)")
            .Define("Tag_polarVec_temp", "ROOT::Math::Polar2DVector(Tag_pt, Tag_phi)")
            .Define("Tag_raw_polarVec_temp", "ROOT::Math::Polar2DVector(Tag_rawPt, Tag_phi)")
            .Define("Probe_polarVec_temp", "ROOT::Math::Polar2DVector(Probe_pt, Probe_phi)")
            .Define("Probe_raw_polarVec_temp", "ROOT::Math::Polar2DVector(Probe_rawPt, Probe_phi)")
            .Define("PuppiMET_polarVec_temp", "ROOT::Math::Polar2DVector(PuppiMET_pt, \
                    PuppiMET_phi)")
            .Define("RawPuppiMET_polarVec_temp", "ROOT::Math::Polar2DVector(RawPuppiMET_pt, \
                    RawPuppiMET_phi)")
    )

    # Activity vector for HDM
    # For dijet and multijet this is the sum of all the jets minus the tag and probe
    # For zjet and egamma this is the sum of all the jets minus the probe
    rdf = (rdf.Define("JetActivity_fourVec_temp", 
                    "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Jet_pt, Jet_eta, \
                            Jet_phi, Jet_mass)")
            .Redefine("JetActivity_fourVec_temp", "ROOT::VecOps::Sum(JetActivity_fourVec_temp, \
                    ROOT::Math::PtEtaPhiMVector())")
    )
    if channel == "dijet" or channel == "multijet":
        rdf = (rdf.Redefine("JetActivity_fourVec_temp",
                            "JetActivity_fourVec_temp - Tag_fourVec_temp - Probe_fourVec_temp"))
    elif channel == "zjet" or channel == "egamma":
        rdf = (rdf.Redefine("JetActivity_fourVec_temp",
                            "JetActivity_fourVec_temp - Probe_fourVec_temp"))
    rdf = (rdf.Define("JetActivity_pt", "float(JetActivity_fourVec_temp.Pt())")
            .Define("JetActivity_eta", "float(JetActivity_fourVec_temp.Eta())")
            .Define("JetActivity_phi", "float(JetActivity_fourVec_temp.Phi())")
            .Define("JetActivity_mass", "float(JetActivity_fourVec_temp.M())")
            .Define("JetActivity_polarVec_temp",
                "ROOT::Math::Polar2DVector(JetActivity_pt, JetActivity_phi)")
    )
    
    return rdf

def do_JEC(rdf):
    rdf = (rdf.Define("DB_direct",
                "-1.0 * Tag_polarVec_temp.Dot(Probe_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("DB_raw_direct",
                "-1.0 * Tag_raw_polarVec_temp.Dot(Probe_raw_polarVec_temp) / (Tag_rawPt * Tag_rawPt)")
            .Define("DB_ratio", "Probe_pt / Tag_pt")
            .Define("DB_raw_ratio", "Probe_rawPt / Tag_rawPt")
            .Define("MPF_tag",
                "1.0 + PuppiMET_polarVec_temp.Dot(Tag_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("MPF_raw_tag",
                "1.0 + RawPuppiMET_polarVec_temp.Dot(Tag_raw_polarVec_temp) / (Tag_rawPt * Tag_rawPt)")
            .Define("MPF_probe",
                "1.0 + PuppiMET_polarVec_temp.Dot(Probe_polarVec_temp) / (Probe_pt * Probe_pt)")
            .Define("MPF_raw_probe",
                "1.0 + RawPuppiMET_polarVec_temp.Dot(Probe_raw_polarVec_temp) / (Probe_rawPt * Probe_rawPt)")
            .Define("R_un_reco_tag_temp",
                "JetActivity_polarVec_temp.Dot(Tag_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("R_un_gen_tag_temp", "1.0")
            .Define("R_un_reco_probe_temp",
                "JetActivity_polarVec_temp.Dot(Probe_polarVec_temp) / (Probe_pt * Probe_pt)")
            .Define("R_un_gen_probe_temp", "1.0")
            .Define("HDM_tag",
                "(DB_direct + MPF_tag - 1.0 + R_un_reco_tag_temp - R_un_gen_tag_temp) / \
                        (cos(ROOT::VecOps::DeltaPhi(Tag_phi, Probe_phi)))")
            .Define("HDM_probe",
                "(DB_direct + MPF_probe - 1.0 + R_un_reco_probe_temp - R_un_gen_probe_temp) / \
                        (cos(ROOT::VecOps::DeltaPhi(Tag_phi, Probe_phi)))")
           )

    # Energy Fraction balance
    rdf = (rdf.Define("EFB_chEmHEF", "(Probe_rawPt * Probe_chEmEF) / Tag_pt")
        .Define("EFB_chHEF", "(Probe_rawPt * Probe_chHEF) / Tag_pt")
        .Define("EFB_hfEmEF", "(Probe_rawPt * Probe_hfEmEF) / Tag_pt")
        .Define("EFB_hfHEF", "(Probe_rawPt * Probe_hfHEF) / Tag_pt")
        .Define("EFB_muEF", "(Probe_rawPt * Probe_muEF) / Tag_pt")
        .Define("EFB_neEmEF", "(Probe_rawPt * Probe_neEmEF) / Tag_pt")
        .Define("EFB_neHEF", "(Probe_rawPt * Probe_neHEF) / Tag_pt")
    )

    return rdf


def run(args):
    # shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    files: List[str] = []
    if args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        for path in paths:
            files.extend(file_read_lines(path))
    else:
        files = [s.strip() for s in args.filelist.split(',')]

    if args.nsteps and args.step:
        n = args.nsteps
        i = args.step
        files = files[i::n]

    triggers: List[str] = []
    if args.triggerlist:
        triggers = args.triggerlist.split(",")
    elif args.triggerpath:
        triggers = file_read_lines(args.triggerpath)

    # Write and hadd the output
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    skim(files, triggers, args, args.step)


def skim(files, triggers, args, step=None):
    # Load the files
    events_chain = ROOT.TChain("Events")
    runs_chain = ROOT.TChain("Runs")

    for file in files:
        if not args.is_local:
            try:
                events_chain.Add(f"root://xrootd-cms.infn.it//{file}")
                runs_chain.Add(f"root://xrootd-cms.infn.it//{file}")
            except Exception as e:
                print(f"Skipping problematic run: {e}")
        else:
            events_chain.Add(file)
            runs_chain.Add(file)

    events_rdf = ROOT.RDataFrame(events_chain)
    runs_rdf = ROOT.RDataFrame(runs_chain)

    if args.progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(events_rdf)

    if args.golden_json:
        events_rdf = filter_json(events_rdf, args.golden_json)

    events_rdf = events_rdf.Filter("nJet > 0", "nJet > 0")

    # Apply corrections
    if args.correction_json:
        import json
        import correctionlib
        correctionlib.register_pyroot_binding()

        with open(args.correction_json) as f:
            correction_info = json.load(f)

        correction_info = correction_info[args.correction_key]

        # JECs
        if "jec_path" and "jec_stack" in correction_info:
            events_rdf = correct_jets(events_rdf, correction_info["jec_path"], correction_info["jec_stack"])
            events_rdf = sort_jets(events_rdf, jet_columns)

        # Vetomaps
        if "vetomap_path" and "vetomap_set" in correction_info:
            events_rdf = find_vetojets(events_rdf, correction_info["vetomap_path"], correction_info["vetomap_set"])
    else:
        # Define that no jets were vetoed
        events_rdf = events_rdf.Define("Jet_vetoed", "ROOT::VecOps::RVec<bool>(Jet_pt.size(), false)")


    # Initialize the JEC variables
    print("Initializing TnP variables")
    events_rdf = init_TnP(events_rdf, args.channel)
    print("Initializing JEC variables")
    events_rdf = do_JEC(events_rdf)


    # Filter based on triggers and one jet
    # Check that the triggers are in the file
    cols = events_rdf.GetColumnNames()
    for trigger in triggers:
        if trigger not in cols and "&&" not in trigger:
            print(f"Trigger {trigger} not in the file") 
            events_rdf = events_rdf.Define(trigger, "0")

    if len(triggers) == 0:
        trg_filter = "1"
    else:
        trg_filter = " || ".join(triggers)
    flag_filter = " && ".join(get_Flags())
    events_rdf = (events_rdf.Filter(trg_filter, trg_filter)
            .Filter(flag_filter, flag_filter)
            )

    # Define a weight column
    if args.is_mc:
        xsec = weight_info["xsec"].get(args.mc_tag)
        if xsec:
            print(f"Reweight with xsec={xsec}")
            events_rdf = (events_rdf.Define("weight", f"{xsec}*genWeight"))
        else:
            events_rdf = (events_rdf.Define("weight", "genWeight"))
    else:
        events_rdf = (events_rdf.Define("weight", "1.0"))

    # Define Pileup_ columns for data
    if not args.is_mc:
        events_rdf = (events_rdf.Define("Pileup_gpudensity", "0.0")
                    .Define("Pileup_nPU", "0")
                    .Define("Pileup_nTrueInt", "0.0")
                    .Define("Pileup_pthatmax", "0.0")
                    .Define("Pileup_sumEOOT", "0.0")
                    .Define("Pileup_sumLOOT", "0.0")
        )

    # Set name of the output file
    step_str = ""
    if step:
        step_str = f"_{step}"
    if args.run_range:
        run_range = args.run_range.split(",")
        assert(len(run_range) == 2)

        print(f"Run range: ({run_range[0]}, {run_range[1]})");
        run_range_str = f"runs{run_range[0]}to{run_range[1]}"

        subprocess.run(["brilcalc", "lumi", "--normtag",
            "/cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json",
            "-u", "/fb", "--begin", f"{run_range[0]}", "--end", f"{run_range[1]}",
            "-i", args.golden_json, "-o", "lumi.csv"])

        df = pd.read_csv("lumi.csv", comment='#', names=["run:fill", "time", "nls",
            "ncms", "delivered(/fb)", "recorded(/fb)"])
        int_lumi = np.sum(df["recorded(/fb)"].to_numpy())
        print(f"Running on {int_lumi} 1/fb integrated luminosity")

        events_rdf = events_rdf.Define("min_run", f"{run_range[0]}")
        events_rdf = events_rdf.Define("max_run", f"{run_range[1]}")
        events_rdf = events_rdf.Define("int_lumi", f"{int_lumi}")
        output_path = os.path.join(args.out, f"J4PSkim_{run_range_str}_{args.channel}{step_str}")
    elif args.mc_tag:
        output_path = os.path.join(args.out, f"J4PSkim_{args.mc_tag}_{args.channel}{step_str}")
        events_rdf = events_rdf.Define("min_run", "0")
        events_rdf = events_rdf.Define("max_run", "1")
        events_rdf = events_rdf.Define("int_lumi", "1.")
    else:
        output_path = os.path.join(args.out, f"J4PSkim_{args.channel}{step_str}")
        events_rdf = events_rdf.Define("min_run", "0")
        events_rdf = events_rdf.Define("max_run", "1")
        events_rdf = events_rdf.Define("int_lumi", "1.")

    # Include defined columns
    columns = [str(col) for col in events_rdf.GetDefinedColumnNames() if not str(col).endswith("_temp") and not str(col).startswith("Jet_")]

    # Include pileup
    columns.extend([str(col) for col in events_rdf.GetColumnNames() if (str(col).startswith("Pileup_") \
                    or str(col).startswith("Rho_") or str(col).startswith("PV_")) and not str(col).endswith("_temp")])
    
    # Include MET
    columns.extend([str(col) for col in events_rdf.GetColumnNames() if (str(col).startswith("RawPFMET") \
                    or str(col).startswith("RawPuppiMET") or str(col).startswith("PuppiMET") or str(col).startswith("PFMET") \
                    or str(col).startswith("CorrT1METJet") or str(col).startswith("RawPFMET")) and not str(col).endswith("_temp")])

    # Include run info
    columns.extend(["weight", "run", "luminosityBlock", "event", "int_lumi", "min_run", "max_run"])

    # Include triggers
    columns.extend([trig for trig in triggers if trig in events_rdf.GetColumnNames()])

    # Check for duplicates
    columns = list(set(columns))
    columns.sort()

    print(f"Writing output for {output_path}.root")
    start = time.time()
    events_rdf.Snapshot("Events", output_path+"_events.root", columns)
    runs_rdf.Snapshot("Runs", output_path+"_runs.root")
    print(f"snapshot finished in {time.time()-start} s for {output_path}.root")

    start = time.time()
    subprocess.run(["hadd", "-f", "-k", output_path+".root",
    output_path+"_events.root", output_path+"_runs.root"])
    os.remove(output_path+"_events.root")
    os.remove(output_path+"_runs.root")
    print(f"hadd finished in {time.time()-start} s for {output_path}.root")

    print(output_path+".root")

    # Get a report of the processing
    report = events_rdf.Report()

    begin = report.begin()
    end = report.end()
    allEntries = 0 if begin == end else begin.__deref__().GetAll()

    # Collect the cuts
    it = begin
    cuts = []
    while it != end:
        ci = it.__deref__()
        cuts.append({ci.GetName(): {"pass": ci.GetPass(), "all": ci.GetAll(), "eff": ci.GetEff(),
                "cumulativeEff": 100.0 * float(ci.GetPass()) / float(allEntries) \
                    if allEntries > 0 else 0.0}})

        it.__preinc__()

    # Create four histograms with alphanumeric bins
    pass_hist = ROOT.TH1D("pass", "pass", len(cuts), 0, len(cuts))
    pass_hist.SetCanExtend(ROOT.TH1.kAllAxes)
    all_hist = ROOT.TH1D("all", "all", len(cuts), 0, len(cuts))
    all_hist.SetCanExtend(ROOT.TH1.kAllAxes)
    eff_hist = ROOT.TH1D("eff", "eff", len(cuts), 0, len(cuts))
    eff_hist.SetCanExtend(ROOT.TH1.kAllAxes)
    cum_eff_hist = ROOT.TH1D("cum_eff", "cum_eff", len(cuts), 0, len(cuts))
    cum_eff_hist.SetCanExtend(ROOT.TH1.kAllAxes)

    for _, cut in enumerate(cuts):
        #print(i, cut)
        for key, value in cut.items():
            pass_hist.Fill(key, value["pass"])
            all_hist.Fill(key, value["all"])
            eff_hist.Fill(key, value["eff"])
            cum_eff_hist.Fill(key, value["cumulativeEff"])

    pass_hist.SetError(np.zeros(len(cuts), dtype=np.float64))
    all_hist.SetError(np.zeros(len(cuts), dtype=np.float64))
    eff_hist.SetError(np.zeros(len(cuts), dtype=np.float64))
    cum_eff_hist.SetError(np.zeros(len(cuts), dtype=np.float64))

    # Save the histograms to test.root
    f = ROOT.TFile(output_path+".root", "UPDATE")
    pass_hist.Write()
    all_hist.Write()
    eff_hist.Write()
    cum_eff_hist.Write()
    f.Close()

