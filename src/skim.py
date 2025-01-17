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

weight_info = {
    "xsec" : {
        # dijet and multijet
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
        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": 6695.0,
        # egamma
        "GJ-4Jets_HT-40to70_TuneCP5_13p6TeV_madgraphMLM-pythia8": 15240.0,
        "GJ-4Jets_HT-70to100_TuneCP5_13p6TeV_madgraphMLM-pythia8": 8111.0,
        "GJ-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8": 7327.0,
        "GJ-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8": 1541.0,
        "GJ-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 167.6,
        "GJ-4Jets_HT-600_TuneCP5_13p6TeV_madgraphMLM-pythia8": 54.39

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
    "Jet_rawFactor"
]

def do_cut_golden_json(rdf, golden_json):
    ROOT.gInterpreter.Declare(
"""
#ifndef GOLDEN_JSON
#define GOLDEN_JSON

#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>

using json = nlohmann::json;

json golden_json;

void init_json(std::string jsonFile) {
    std::cout << "Initializing JSON file" << std::endl;
    std::ifstream f(jsonFile);
    golden_json = json::parse(f);
}

bool isGoodLumi(int run, int lumi) {
   for (auto& lumiRange : golden_json[std::to_string(run)]) {
       if (lumi >= lumiRange[0] && lumi <= lumiRange[1]) {
           return true;
       }
   }

    return false;
}

#endif
"""
    )
    ROOT.init_json(golden_json)
    print("Applying golden JSON cut")
    print(f"JSON file: {golden_json}")
    rdf = (rdf.Filter("isGoodLumi(run, luminosityBlock)", "Golden JSON"))
    return rdf

def init_TnP(rdf, dataset):
    # TODO:
    # - Implement an index finding function for probe jets, to increase clarity
    #   and to avoid cutting on the jet collection
    if dataset == "dijet":
        ROOT.gInterpreter.Declare("""
        #ifndef DIJET_IDXS
        #define DIJET_IDXS
                                  
        std::pair<int, int> findTagProbeIdxs(ROOT::RVec<float> Jet_eta, ROOT::RVec<float> Jet_pt,
                            ROOT::RVec<float> Jet_phi, ROOT::RVec<int> Jet_jetId) {
            int idx1 = -1;
            int idx2 = -1;
                                  
            // Find the tag jet as the leading barrel jet
            for (int i = 0; i < Jet_pt.size(); i++) {
                if (abs(Jet_eta[i]) < 1.3 && Jet_pt[i] > 15 && Jet_jetId[i] >= 4) {
                    idx1 = i;
                    break;
                }
            }

            // Find the probe jet as:
            // leading jet back-to-back with the tag jet
            for (int i = 0; i < Jet_pt.size(); i++) {
                if (i == idx1 || Jet_pt[i] < 15) {
                    continue;
                }
                if (abs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Jet_phi[idx1])) > 2.7 &&
                    Jet_pt[i] > 12 && Jet_jetId[i] >= 4) {
                    idx2 = i;
                    break;
                }
            }

            // If probe also in barrel, randomize the indices
            if (idx2 != -1 && abs(Jet_eta[idx2]) < 1.3 && idx1 != -1) {
                bool swap = (int(Jet_phi[idx1] * 100) % 2) == 0; // Jets are ~uniform in phi, \
                        use this to generate a random number
                if (swap) {
                    std::swap(idx1, idx2);
                }
            }

            return std::make_pair(idx1, idx2);
        }
                                  
        #endif
        """)

        rdf = (rdf.Define("TnP_idx_temp", "findTagProbeIdxs(Jet_eta, Jet_pt, Jet_phi, Jet_jetId)")
                .Filter("TnP_idx_temp.first >= 0 && TnP_idx_temp.second >= 0",
                    "Tag and probe found")
                .Define("Tag_idx_temp", "TnP_idx_temp.first")
                .Define("Probe_idx_temp", "TnP_idx_temp.second")
                .Define("Tag_pt", "Jet_pt[Tag_idx_temp]")
                .Define("Tag_eta", "Jet_eta[Tag_idx_temp]")
                .Define("Tag_phi", "Jet_phi[Tag_idx_temp]")
                .Define("Tag_mass", "Jet_mass[Tag_idx_temp]")
                .Define("Tag_label", "0")
        )

        # Create a probe jet collection
        for column in jet_columns:
            rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")


    elif dataset == "zjet":
        ROOT.gInterpreter.Declare("""
        #ifndef ZJET_IDXS
        #define ZJET_IDXS
                                  
        std::pair<int, int> findMuonIdxs(ROOT::RVec<float> Muon_eta, ROOT::RVec<float> Muon_pt,
                            ROOT::RVec<float> Muon_pfRelIso03_all, ROOT::RVec<int> Muon_tightId,
                            ROOT::RVec<float> Muon_charge) {
            int idx1 = -1;
            int idx2 = -1;
            for (int i = 0; i < Muon_pt.size(); i++) {
                if (abs(Muon_eta[i]) < 2.4 && Muon_pfRelIso03_all[i] < 0.15 &&
                    Muon_tightId[i]) {
                    // Leading muon pt>20, subleading pt>10
                    if (idx1 == -1 && Muon_pt[i] > 20) {
                        idx1 = i;
                    } else if (idx2 == -1 && Muon_charge[i] != Muon_charge[idx1] &&
                            Muon_pt[i] > 10) {
                        idx2 = i;
                        break;
                    }
                }
            }
                                  
            return std::make_pair(idx1, idx2);
        }

        int findJetIdx(ROOT::RVec<float> Jet_eta, ROOT::RVec<float> Jet_pt,
                            ROOT::RVec<float> Jet_phi, ROOT::RVec<int> Jet_jetId,
                            float Z_eta, float Z_phi) {
            for (int i = 0; i < Jet_pt.size(); i++) {
                if (Jet_pt[i] > 12 && Jet_jetId[i] >= 4 && abs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Z_phi)) > 2.7) {
                    return i;
                }
            }
                                  
            return -1;
        }
                                  
        #endif
        """)
        rdf = (rdf.Filter("nMuon > 1", "nMuon > 1")
                .Define("Muon_idx_temp", "findMuonIdxs(Muon_eta, Muon_pt, Muon_pfRelIso03_all, \
                        Muon_tightId, Muon_charge)")
                .Filter("Muon_idx_temp.first >= 0 && Muon_idx_temp.second >= 0", "Two muons found")
                .Define("Z_4vec_temp",
                    "ROOT::Math::PtEtaPhiMVector(Muon_pt[Muon_idx_temp.first], \
                            Muon_eta[Muon_idx_temp.first], Muon_phi[Muon_idx_temp.first], \
                            Muon_mass[Muon_idx_temp.first]) + \
                            ROOT::Math::PtEtaPhiMVector(Muon_pt[Muon_idx_temp.second], \
                            Muon_eta[Muon_idx_temp.second], Muon_phi[Muon_idx_temp.second], \
                            Muon_mass[Muon_idx_temp.second])")
                .Define("Tag_pt", "static_cast<float>(Z_4vec_temp.Pt())")
                .Define("Tag_eta", "static_cast<float>(Z_4vec_temp.Eta())")
                .Define("Tag_phi", "static_cast<float>(Z_4vec_temp.Phi())")
                .Define("Tag_mass", "static_cast<float>(Z_4vec_temp.M())")
                .Define("Tag_label", "1")
                .Define("Probe_idx_temp", "findJetIdx(Jet_eta, Jet_pt, Jet_phi, Jet_jetId, \
                        Tag_eta, Tag_phi)")
                .Filter("Probe_idx_temp >= 0", "Jet found")
                .Filter("Tag_pt > 12", "Z pT > 12")
                .Filter("Tag_mass > 71.1876 && Tag_mass < 111.1876", "Z mass window")
        )

        for column in jet_columns:
            rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    elif dataset == "egamma":
        ROOT.gInterpreter.Declare("""
        #ifndef EGAMMA_IDXS
        #define EGAMMA_IDXS
        
        int findPhotonIdx(ROOT::RVec<float> Photon_eta, ROOT::RVec<float> Photon_pt,
                            ROOT::RVec<int> Photon_cutBased, ROOT::RVec<float> Photon_hoe,
                            ROOT::RVec<float> Photon_r9) {
            for (int i = 0; i < Photon_pt.size(); i++) {
                if (abs(Photon_eta[i]) < 1.3 && Photon_pt[i] > 15 && Photon_cutBased[i] == 3 &&
                    Photon_hoe[i] < 0.02148 && Photon_r9[i] > 0.94 && Photon_r9[i] < 1.00) {
                    return i;
                }
            }
                                  
            return -1;
        }
                                  
        int findJetIdx(ROOT::RVec<float> Jet_eta, ROOT::RVec<float> Jet_pt,
                            ROOT::RVec<float> Jet_phi, ROOT::RVec<int> Jet_jetId,
                            int Photon_jetIdx, float Photon_phi) {
            for (int i = 0; i < Jet_pt.size(); i++) {
                if (Jet_pt[i] > 12 && Jet_jetId[i] >= 4
                    && i != Photon_jetIdx &&
                    abs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Photon_phi)) > 2.7) {
                    return i;
                }
            }
                                  
            return -1;
        }
                                  
        #endif
        """)

        rdf = (rdf.Define("Tag_idx_temp", "findPhotonIdx(Photon_eta, Photon_pt, Photon_cutBased, \
                Photon_hoe, Photon_r9)")
                .Filter("Tag_idx_temp >= 0", "Photon found")
                .Define("Probe_idx_temp", "findJetIdx(Jet_eta, Jet_pt, Jet_phi, Jet_jetId, \
                        Photon_jetIdx[Tag_idx_temp], Photon_phi[Tag_idx_temp])")
                .Filter("Probe_idx_temp >= 0", "Jet found")
                .Define("Tag_pt", f"Photon_pt[Tag_idx_temp]")
                .Define("Tag_eta", f"Photon_eta[Tag_idx_temp]")
                .Define("Tag_phi", f"Photon_phi[Tag_idx_temp]")
                .Define("Tag_mass", "0.0")
                .Define("Tag_label", "2")
        )

        for column in jet_columns:
            rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    elif dataset == "multijet":
        # Change Tag <-> Probe for multijet, since low pt jets better calibrated?
        recoil_filter = "abs(RecoilJet_eta)<2.5 && RecoilJet_pt>30"
        rdf = (rdf.Filter("nJet > 2", "nJet > 2")
                .Filter("Jet_pt[0] > 30",
                    "Leading jet pT > 30")
                .Define("Probe_pt", "Jet_pt[0]")
                .Define("Probe_eta", "Jet_eta[0]")
                .Define("Probe_phi", "Jet_phi[0]")
                .Define("Probe_mass", "Jet_mass[0]")
                .Define("Tag_label", "3")
                .Define("RecoilJet_ids",
                    "ROOT::VecOps::Drop(ROOT::VecOps::Enumerate(Jet_pt), \
                            ROOT::VecOps::RVec<int>{0})")
                .Define("RecoilJet_pt", f"ROOT::VecOps::Take(Jet_pt, RecoilJet_ids)")
                .Define("RecoilJet_eta", f"ROOT::VecOps::Take(Jet_eta, RecoilJet_ids)")
                .Define("RecoilJet_phi", f"ROOT::VecOps::Take(Jet_phi, RecoilJet_ids)")
                .Define("RecoilJet_mass", f"ROOT::VecOps::Take(Jet_mass, RecoilJet_ids)")
                .Define("Tag_pt", f"RecoilJet_pt[{recoil_filter}]")
                .Define("Tag_eta", f"RecoilJet_eta[{recoil_filter}]")
                .Define("Tag_phi", f"RecoilJet_phi[{recoil_filter}]")
                .Define("Tag_mass", f"RecoilJet_mass[{recoil_filter}]")
                .Define("Tag_ids", f"RecoilJet_ids[{recoil_filter}]")
                .Define("TagMJ_fourVec_temp",
                    "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Tag_pt, \
                            Tag_eta, Tag_phi, Tag_mass)")
                .Redefine("TagMJ_fourVec_temp",
                    "ROOT::VecOps::Sum(TagMJ_fourVec_temp, ROOT::Math::PtEtaPhiMVector())")
                .Redefine("Tag_pt",
                    "float(TagMJ_fourVec_temp.Pt())")
                .Redefine("Tag_eta",
                    "float(TagMJ_fourVec_temp.Eta())")
                .Redefine("Tag_phi",
                    "float(TagMJ_fourVec_temp.Phi())")
                .Redefine("Tag_mass",
                    "float(TagMJ_fourVec_temp.M())")
        )

        for column in jet_columns:
            if column in ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"]:
                continue
            # For multijet change Probe columns to be zero, as probe is not a jet
            rdf = rdf.Define("Probe_"+column[4:], f"0.0")

    rdf = rdf.Define("Probe_rawPt", "(1.0 - Probe_rawFactor) * Probe_pt")

    # Label non-flat branches as _temp to drop them later
    rdf = (rdf.Define("Tag_fourVec_temp", "ROOT::Math::PtEtaPhiMVector(Tag_pt, Tag_eta, Tag_phi, Tag_mass)")
            .Define("Probe_fourVec_temp", "ROOT::Math::PtEtaPhiMVector(Probe_pt, Probe_eta, \
                    Probe_phi, Probe_mass)")
            .Define("Tag_polarVec_temp", "ROOT::Math::Polar2DVector(Tag_pt, Tag_phi)")
            .Define("Probe_polarVec_temp", "ROOT::Math::Polar2DVector(Probe_pt, Probe_phi)")
            .Define("PuppiMET_polarVec_temp", "ROOT::Math::Polar2DVector(PuppiMET_pt, \
                    PuppiMET_phi)")
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
    if dataset == "dijet" or dataset == "multijet":
        rdf = (rdf.Redefine("JetActivity_fourVec_temp",
                            "JetActivity_fourVec_temp - Tag_fourVec_temp - Probe_fourVec_temp"))
    elif dataset == "zjet" or dataset == "egamma":
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
            .Define("DB_ratio", "Probe_pt / Tag_pt")
            .Define("MPF_tag",
                "1.0 + PuppiMET_polarVec_temp.Dot(Tag_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("MPF_probe",
                "1.0 + PuppiMET_polarVec_temp.Dot(Probe_polarVec_temp) / (Probe_pt * Probe_pt)")
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

    files: List[str] = []
    if args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        for path in paths:
            files.extend(file_read_lines(path))
    else:
        files = [s.strip() for s in args.filelist.split(',')]

    triggers: List[str] = []
    if args.triggerlist:
        triggers = args.triggerlist.split(",")
    elif args.triggerpath:
        triggers = file_read_lines(args.triggerpath)

    # Write and hadd the output
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    # Load the files
    print(f"Processing files")
    events_chain = ROOT.TChain("Events")
    runs_chain = ROOT.TChain("Runs")

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

    if args.golden_json:
        events_rdf = do_cut_golden_json(events_rdf, args.golden_json)

    events_rdf = events_rdf.Filter("nJet > 0", "nJet > 0")

    # Initialize the JEC variables
    print("Initializing TnP variables")
    events_rdf = init_TnP(events_rdf, args.dataset)
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

    # Set name of the output file
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
        output_path = os.path.join(args.out, f"J4PSkim_{run_range_str}_{args.dataset}")
    elif args.mc_tag:
        output_path = os.path.join(args.out, f"J4PSkim_{args.mc_tag}_{args.dataset}")
        events_rdf = events_rdf.Define("min_run", "0")
        events_rdf = events_rdf.Define("max_run", "1")
        events_rdf = events_rdf.Define("int_lumi", "1.")
    else:
        output_path = os.path.join(args.out, f"J4PSkim_{args.dataset}")
        events_rdf = events_rdf.Define("min_run", "0")
        events_rdf = events_rdf.Define("max_run", "1")
        events_rdf = events_rdf.Define("int_lumi", "1.")


    # Remove the Jet_ and _temp columns
    if args.defined_columns:
        pass


    # Keep only the columns that are needed
    columns = [str(col) for col in events_rdf.GetColumnNames() if (str(col).startswith("Probe_") or str(col).startswith("Tag_") \
               or "DB_" in str(col) or "MPF_" in str(col) or "HDM_" in str(col) or "EFB_" in str(col)) and not str(col).endswith("_temp")]

    columns.extend(["weight", "run", "luminosityBlock", "event", "int_lumi", "min_run", "max_run"])
    columns.extend([trig for trig in triggers if trig in events_rdf.GetColumnNames()])


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

