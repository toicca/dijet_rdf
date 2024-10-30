import ROOT
import os
import subprocess
import argparse
import pathlib
import ctypes
import numpy as np

from RDFHelpers import file_read_lines
from processing_utils import find_site

jet_columns = [
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId",
    "Jet_area", "Jet_nConstituents", "Jet_nElectrons", "Jet_nMuons",
    "Jet_chEmEF", "Jet_neEmEF", "Jet_chHEF", "Jet_neHEF",
    "Jet_rawFactor"
]

def init_TnP(rdf, dataset):
    if dataset == "dijet":
        rdf = (rdf.Filter("nJet > 1", "nJet > 1")
                .Filter("abs(ROOT::VecOps::DeltaPhi(Jet_phi[0], Jet_phi[1])) > 2.7", \
                        "abs(DeltaPhi(Jet1, Jet2)) > 2.7")
                .Filter("nJet > 2 ? Jet_pt[2] / ((Jet_pt[0] + Jet_pt[1]) * 0.5) < 1.0 : true", "alpha < 1.0")
                .Filter("(Jet_pt[0]/Jet_pt[1] < 1.3 && Jet_pt[0]/Jet_pt[1] > 0.7)", "1.3 > pT1/pT2 > 0.7")
                .Define("firstTag", "int(abs(Jet_eta[0]) < 1.3)")
                .Define("secondTag", "int(abs(Jet_eta[1]) < 1.3)")
                .Define("nTag", "firstTag + secondTag")
                .Filter("nTag > 0", "nTag > 0")
                .Define("Tag_ids", "nTag == 2 ? ROOT::VecOps::RVec<int>{0, 1} : \
                        (firstTag == 1 ? ROOT::VecOps::RVec<int>{0} : \
                        ROOT::VecOps::RVec<int>{1})")
                .Define("Probe_ids", "1 - Tag_ids")
                .Define("Tag_pt", "ROOT::VecOps::Take(Jet_pt, Tag_ids)")
                .Define("Tag_eta", "ROOT::VecOps::Take(Jet_eta, Tag_ids)")
                .Define("Tag_phi", "ROOT::VecOps::Take(Jet_phi, Tag_ids)")
                .Define("Tag_mass", "ROOT::VecOps::Take(Jet_mass, Tag_ids)")
                .Define("Tag_label", "ROOT::VecOps::RVec<int>(Tag_pt.size(), 0)")
                .Define("JetT_ids_temp", "ROOT::VecOps::Enumerate(Jet_pt)")
                .Define("TnP_ids_temp", "ROOT::VecOps::Concatenate(Tag_ids, Probe_ids)")
                .Define("JetActivity_ids", "ROOT::VecOps::Drop(JetT_ids_temp, TnP_ids_temp)")
        )
        for column in jet_columns:
            rdf = rdf.Define("Probe_"+column[4:], f"ROOT::VecOps::Take({column}, Probe_ids)")

    elif dataset == "zjet":
        muon_filter = "abs(Muon_eta)<2.4 && Muon_pt>20 && Muon_pfRelIso04_all<0.15 && Muon_tightId"
        rdf = (rdf.Filter("nJet > 0", "nJet > 0")
                .Filter("nMuon > 1", "nMuon > 1")
                .Filter("Muon_charge[0] + Muon_charge[1] == 0", "Opposite charge muons")
                .Define("ZMuons_pt", f"Muon_pt[{muon_filter}]")
                .Define("ZMuons_eta", f"Muon_eta[{muon_filter}]")
                .Define("ZMuons_phi", f"Muon_phi[{muon_filter}]")
                .Define("ZMuons_mass", f"Muon_mass[{muon_filter}]")
                .Filter("ZMuons_pt.size() == 2", "Exactly 2 muons")
                .Define("Z_4vec",
                    "ROOT::Math::PtEtaPhiMVector(ZMuons_pt[0], ZMuons_eta[0], ZMuons_phi[0], \
                            ZMuons_mass[0]) + ROOT::Math::PtEtaPhiMVector(ZMuons_pt[1], \
                            ZMuons_eta[1], ZMuons_phi[1], ZMuons_mass[1])")
                .Define("Tag_pt", "ROOT::VecOps::RVec<float>{static_cast<float>(Z_4vec.Pt())}")
                .Define("Tag_eta", "ROOT::VecOps::RVec<float>{static_cast<float>(Z_4vec.Eta())}")
                .Define("Tag_phi", "ROOT::VecOps::RVec<float>{static_cast<float>(Z_4vec.Phi())}")
                .Define("Tag_mass", "ROOT::VecOps::RVec<float>{static_cast<float>(Z_4vec.M())}")
                .Define("Tag_label", "ROOT::VecOps::RVec<int>{1}")
                .Define("nTag", "Tag_pt.size()")
                .Define("Probe_ids", "ROOT::VecOps::RVec<int>{0}")
                .Define("JetActivity_ids",
                    "ROOT::VecOps::Drop(ROOT::VecOps::Enumerate(Jet_pt), Probe_ids)")
                .Filter("abs(ROOT::VecOps::DeltaPhi(Jet_phi[0], Tag_phi[0])) > 2.7",
                    "abs(DeltaPhi(Jet1, Z)) > 2.7")
                .Filter("nJet > 1 ? Jet_pt[1] / Tag_pt[0] < 1.0 : true", "alpha < 1.0")
        )

        for column in jet_columns:
            rdf = rdf.Define("Probe_"+column[4:], f"ROOT::VecOps::Take({column}, 1)")

    elif dataset == "egamma":
        photon_filter = "abs(Photon_eta)<1.3 && Photon_pt>15"
        rdf = (rdf.Filter("nJet > 0", "nJet > 0")
                .Filter("nPhoton > 0", "nPhoton > 0")
                .Define("Tag_pt", f"ROOT::VecOps::Take(Photon_pt[{photon_filter}], 1, (float)0.0)")
                .Define("Tag_eta", f"ROOT::VecOps::Take(Photon_eta[{photon_filter}], 1, (float)0.0)")
                .Define("Tag_phi", f"ROOT::VecOps::Take(Photon_phi[{photon_filter}], 1, (float)0.0)")
                .Define("Tag_mass", "ROOT::VecOps::RVec<float>(1, 0.0)")
                .Define("Tag_label", "ROOT::VecOps::RVec<int>{2}")
                .Define("nTag", "Tag_pt.size()")
                .Define("Probe_ids", "ROOT::VecOps::RVec<int>{0}")
                .Define("JetActivity_ids", "ROOT::VecOps::Drop(ROOT::VecOps::Enumerate(Jet_pt), Probe_ids)")
                .Filter("nTag == 1", "Exactly 1 photon")
                .Filter("abs(ROOT::VecOps::DeltaPhi(Jet_phi[0], Tag_phi[0])) > 2.7",
                        "abs(DeltaPhi(Jet1, Photon)) > 2.7")
                .Filter("nJet > 1 ? Jet_pt[1] / Tag_pt[0] < 1.0 : true", "alpha < 1.0")
        )

        for column in jet_columns:
            rdf = rdf.Define("Probe_"+column[4:], f"ROOT::VecOps::Take({column}, 1)")

    elif dataset == "multijet":
        # Change Tag <-> Probe for multijet, since low pt jets better calibrated?
        recoil_filter = "abs(RecoilJet_eta)<2.5 && RecoilJet_pt>30"
        rdf = (rdf.Filter("nJet > 2", "nJet > 2")
                .Filter("Jet_pt[0] > 30 && abs(Jet_eta[0]) < 2.5", "Leading jet")
                .Define("Tag_pt", "ROOT::VecOps::RVec<float>{Jet_pt[0]}")
                .Define("Tag_eta", "ROOT::VecOps::RVec<float>{Jet_eta[0]}")
                .Define("Tag_phi", "ROOT::VecOps::RVec<float>{Jet_phi[0]}")
                .Define("Tag_mass", "ROOT::VecOps::RVec<float>{Jet_mass[0]}")
                .Define("Tag_label", "ROOT::VecOps::RVec<int>{3}")
                .Define("nTag", "Tag_pt.size()")
                .Define("RecoilJet_ids",
                    "ROOT::VecOps::Drop(ROOT::VecOps::Enumerate(Jet_pt), \
                            ROOT::VecOps::RVec<int>{0})")
                .Define("RecoilJet_pt", f"ROOT::VecOps::Take(Jet_pt, RecoilJet_ids)")
                .Define("RecoilJet_eta", f"ROOT::VecOps::Take(Jet_eta, RecoilJet_ids)")
                .Define("RecoilJet_phi", f"ROOT::VecOps::Take(Jet_phi, RecoilJet_ids)")
                .Define("RecoilJet_mass", f"ROOT::VecOps::Take(Jet_mass, RecoilJet_ids)")
                .Define("Probe_pt", f"RecoilJet_pt[{recoil_filter}]")
                .Define("Probe_eta", f"RecoilJet_eta[{recoil_filter}]")
                .Define("Probe_phi", f"RecoilJet_phi[{recoil_filter}]")
                .Define("Probe_mass", f"RecoilJet_mass[{recoil_filter}]")
                .Define("Probe_ids", f"RecoilJet_ids[{recoil_filter}]")
                .Define("ProbeMJ_fourVec_temp",
                    "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Probe_pt, \
                            Probe_eta, Probe_phi, Probe_mass)")
                .Redefine("ProbeMJ_fourVec_temp",
                    "ROOT::VecOps::Sum(ProbeMJ_fourVec_temp, ROOT::Math::PtEtaPhiMVector())")
                .Redefine("Probe_pt",
                    "ROOT::VecOps::RVec<float>{float(ProbeMJ_fourVec_temp.Pt())}")
                .Redefine("Probe_eta",
                    "ROOT::VecOps::RVec<float>{float(ProbeMJ_fourVec_temp.Eta())}")
                .Redefine("Probe_phi",
                    "ROOT::VecOps::RVec<float>{float(ProbeMJ_fourVec_temp.Phi())}")
                .Redefine("Probe_mass",
                    "ROOT::VecOps::RVec<float>{float(ProbeMJ_fourVec_temp.M())}")
                .Define("JetActivity_ids", "ROOT::VecOps::Drop(RecoilJet_ids, Probe_ids)")
        )

        for column in jet_columns:
            if column in ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"]:
                continue
            # For multijet change Probe columns to be zero, as probe is not a jet
            rdf = rdf.Define("Probe_"+column[4:], f"ROOT::VecOps::RVec<float>(1, 0.0)")

    # Label non-flat branches as _temp to drop them later
    rdf = (rdf.Define("nProbe", "Probe_pt.size()")
            .Filter("nProbe > 0", "nProbe > 0")
            .Define("Tag_fourVec_temp",
                "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Tag_pt, Tag_eta, \
                        Tag_phi, Tag_mass)")
            .Define("Probe_fourVec_temp",
                "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Probe_pt, Probe_eta, \
                        Probe_phi, Probe_mass)")
            .Define("Tag_polarVec_temp", 
                "ROOT::VecOps::Construct<ROOT::Math::Polar2DVector>(Tag_pt, Tag_phi)")
            .Define("Probe_polarVec_temp",
                "ROOT::VecOps::Construct<ROOT::Math::Polar2DVector>(Probe_pt, Probe_phi)")
            .Define("PuppiMET_polarVec_probe_temp",
                "ROOT::VecOps::RVec<ROOT::Math::Polar2DVector>(nProbe, \
                        ROOT::Math::Polar2DVector(PuppiMET_pt, PuppiMET_phi))")
            .Define("PuppiMET_polarVec_tag_temp",
                "ROOT::VecOps::RVec<ROOT::Math::Polar2DVector>(nTag, \
                        ROOT::Math::Polar2DVector(PuppiMET_pt, PuppiMET_phi))")
    )

    # Activity vector for HDM
    rdf = (rdf.Define("JetActivity_pt", "ROOT::VecOps::Take(Jet_pt, JetActivity_ids)")
            .Define("JetActivity_eta", "ROOT::VecOps::Take(Jet_eta, JetActivity_ids)")
            .Define("JetActivity_phi", "ROOT::VecOps::Take(Jet_phi, JetActivity_ids)")
            .Define("JetActivity_mass", "ROOT::VecOps::Take(Jet_mass, JetActivity_ids)")
            .Define("JetActivity_fourVec_temp", 
                "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(JetActivity_pt, \
                        JetActivity_eta, JetActivity_phi, JetActivity_mass)")
            .Redefine("JetActivity_fourVec_temp", "ROOT::VecOps::Sum(JetActivity_fourVec_temp, \
                    ROOT::Math::PtEtaPhiMVector())")
            .Redefine("JetActivity_pt", "float(JetActivity_fourVec_temp.Pt())")
            .Redefine("JetActivity_eta", "float(JetActivity_fourVec_temp.Eta())")
            .Redefine("JetActivity_phi", "float(JetActivity_fourVec_temp.Phi())")
            .Redefine("JetActivity_mass", "float(JetActivity_fourVec_temp.M())")
            .Define("JetActivity_polarVec_tag_temp",
                "ROOT::VecOps::RVec<ROOT::Math::Polar2DVector>(nTag, \
                        ROOT::Math::Polar2DVector(JetActivity_pt, JetActivity_phi))")
            .Define("JetActivity_polarVec_probe_temp", 
                "ROOT::VecOps::RVec<ROOT::Math::Polar2DVector>(nProbe, \
                        ROOT::Math::Polar2DVector(JetActivity_pt, JetActivity_phi))")
    )
    
    return rdf

def do_JEC(rdf):
    ROOT.gInterpreter.Declare('''
#ifndef __JECDOT__
#define __JECDOT__
auto vec2D_dot = [](ROOT::Math::Polar2DVector v1, ROOT::Math::Polar2DVector v2) \
        { return v1.Dot(v2); };
#endif
''')
    rdf = (rdf.Define("DB_direct", "-1.0 * ROOT::VecOps::Map(Tag_polarVec_temp, \
            Probe_polarVec_temp, vec2D_dot) / (Tag_pt * Tag_pt)")
           .Define("DB_ratio", "Probe_pt / Tag_pt")
           .Define("MPF_tag", "1.0 + ROOT::VecOps::Map(PuppiMET_polarVec_tag_temp, \
                   Tag_polarVec_temp, vec2D_dot) / (Tag_pt * Tag_pt)")
           .Define("MPF_probe", "1.0 + ROOT::VecOps::Map(PuppiMET_polarVec_probe_temp, \
                   Probe_polarVec_temp, vec2D_dot) / (Probe_pt * Probe_pt)")
           .Define("HDM_tag", "MPF_tag - ROOT::VecOps::Map(Tag_polarVec_temp, \
                   JetActivity_polarVec_tag_temp, vec2D_dot) / (Tag_pt * Tag_pt)")
           .Define("HDM_probe", "MPF_probe - ROOT::VecOps::Map(Probe_polarVec_temp, \
                   JetActivity_polarVec_tag_temp, vec2D_dot) / (Probe_pt * Probe_pt)")
           )

    return rdf

def get_Flags(campaign=None):
    # TODO: Implement campaign-specific flags
    flags = [
            "Flag_goodVertices",
            "Flag_globalSuperTightHalo2016Filter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter",
            "Flag_BadPFMuonDzFilter",
            "Flag_hfNoisyHitsFilter",
            "Flag_eeBadScFilter",
            "Flag_ecalBadCalibFilter"
    ]

    return flags

def run(args):
    # shut up ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    if args.nThreads:
        ROOT.EnableImplicitMT(args.nThreads)

    events_chain = ROOT.TChain("Events")
    runs_chain = ROOT.TChain("Runs")

    files: List[str] = []
    if args.filepath:
        files = file_read_lines(args.filepath)
    else:
        files = [s.strip() for s in args.filelist.split(',')]
    
    triggers: List[str] = []
    if args.triggerlist:
        triggers = args.triggerlist.split(",")
    elif args.triggerpath:
        triggers = file_read_lines(args.triggerpath)

    run_range = args.run_range.split(",")
    assert(len(run_range) == 2)

    # Load the files
    print("Loading files")

    all_columns = []
    for file in files:
        if not args.is_local:
            events_chain.Add(f"root://cms-xrd-global.cern.ch/{file}")
            runs_chain.Add(f"root://cms-xrd-global.cern.ch/{file}")
        else:
            events_chain.Add(file)
            runs_chain.Add(file)

    events_rdf = ROOT.RDataFrame(events_chain)
    runs_rdf = ROOT.RDataFrame(runs_chain)

    if args.progress_bar:
        ROOT.RDF.Experimental.AddProgressBar(events_rdf)

    # Filter based on triggers and one jet
    trg_filter = " || ".join(triggers)
    flag_filter = " && ".join(get_Flags())
    events_rdf = (events_rdf.Filter(trg_filter, trg_filter)
                .Filter(flag_filter, flag_filter)
                .Filter("nJet > 0", "nJet > 0")
                )

    # Keep tight jetId jets
    for col in jet_columns:
        if not str(col).startswith("Jet_"):
            continue
        if "Jet_jetId" in str(col):
            continue
        events_rdf = events_rdf.Redefine(f"{col}", f"{col}[Jet_jetId == 6]")

    events_rdf = (events_rdf.Redefine("Jet_jetId", "Jet_jetId[Jet_jetId == 6]")
                            .Redefine("nJet", "Jet_pt.size()"))

    # Initialize the JEC variables
    print("Initializing TnP variables")
    events_rdf = init_TnP(events_rdf, args.dataset)
    print("Initializing JEC variables")
    events_rdf = do_JEC(events_rdf)

    # Define a weight column
    events_rdf = events_rdf.Define("weight", "1.0")

    # Remove the Jet_ and _temp columns
    print("Removing unnecessary columns")
    all_columns = events_rdf.GetColumnNames()
    #all_columns.extend(events_rdf.GetDefinedColumnNames())
    all_columns = [str(col) for col in all_columns \
                    if not str(col).startswith("Jet_") and not str(col).endswith("_temp")]

    # Write and hadd the output
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    snapshot_opts = ROOT.RDF.RSnapshotOptions()
    snapshot_opts.fCompressionLevel = 9
    print(f"Run range: ({run_range[0]}, {run_range[1]})");
    output_path = os.path.join(args.out,
            f"J4PSkim_runs{run_range[0]}to{run_range[1]}_{args.run_tag}")

    print("Writing output")
    events_rdf.Snapshot("Events", output_path+"_events.root", all_columns)
    runs_rdf.Snapshot("Runs", output_path+"_runs.root")

    subprocess.run(["hadd", "-f7", output_path+".root",
        output_path+"_events.root", output_path+"_runs.root"])
    os.remove(output_path+"_events.root")
    os.remove(output_path+"_runs.root")

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

    for i, cut in enumerate(cuts):
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

