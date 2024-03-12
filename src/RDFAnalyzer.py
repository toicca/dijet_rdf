import ROOT
from dataclasses import dataclass
from typing import List
from RDFHelpers import get_bins
import numpy as np
from make_JEC import compile_JEC, load_JEC, clean_JEC
    
RDataFrame = ROOT.RDataFrame
RNode = ROOT.RDF.RNode
RunGraphs = ROOT.RDF.RunGraphs

@dataclass
class JEC_corrections:
    L1 : str = ""
    L2Relative : str = ""
    L2L3 : str = ""
    JER : str = ""
    JERSF : str = ""
    
    def check_empty(self):
        return all(value == "" or value is None for value in self.__dict__.values())
    
class RDFAnalyzer:
    def __init__(self, filelist : List[str],
                trigger_list : List[str],
                json_file : str = "",
                nFiles : int = -1,
                JEC : JEC_corrections = JEC_corrections("", "", ""),
                nThreads : int = 1,
                progress_bar : bool = False,
                isMC : bool = False,
                ) -> "RDFAnalyzer":
        self.nThreads = nThreads
        self.trigger_list = trigger_list
        self.histograms = {"all" : []} # format : {trigger : [histograms]}
        self.trigger_rdfs = {} # format : {trigger : rdf}. self.rdf is not initialized here due to order of operations
        self.has_run = False # possibly unnecessary
        self.chain = None
        self.bins = get_bins()
        
        self.rdf = self.__loadRDF(filelist, nFiles = nFiles)
        if progress_bar:
            ROOT.RDF.Experimental.AddProgressBar(self.rdf)

        # Initial variables
        self.rdf = (self.rdf.Define("weight", "genWeight" if isMC else "1.0")
                    .Define("Jet_order", "ROOT::VecOps::Argsort(Jet_pt)")
                    .Define("Jet_rawPt", "Jet_pt * (1.0-Jet_rawFactor)")
                    .Define("Jet_passesVetomap", "ROOT::VecOps::RVec<int>(Jet_pt.size(), 1)")
                    .Define("RawPuppiMET_polar", "ROOT::Math::Polar2DVectorF(RawPuppiMET_pt, RawPuppiMET_phi)")
                    )
        
        if not JEC.check_empty():
            print(JEC.L1, JEC.L2Relative, JEC.L2L3)
            # clean_JEC()
            # compile_JEC()
            load_JEC()
            self.rdf = self.__redo_JEC(JEC)

        if json_file != "":
            self.rdf = self.__do_cut_golden_json(json_file)

        for trigger in trigger_list:
            # Consider changing these to a static size or Pythons array
            self.trigger_rdfs[trigger] = self.rdf.Filter(trigger)
            self.histograms[trigger] = []

        # Only after JECs, common cuts etc. have been applied, we can create the inclusive rdf
        self.trigger_rdfs["all"] = self.rdf
        
        return self

    def Flag_cut(self, rdf : RNode) -> RNode:
        flag = """Flag_goodVertices && 
                    Flag_globalSuperTightHalo2016Filter &&
                    Flag_EcalDeadCellTriggerPrimitiveFilter &&
                    Flag_BadPFMuonFilter &&
                    Flag_BadPFMuonDzFilter && 
                    Flag_hfNoisyHitsFilter &&
                    Flag_eeBadScFilter &&
                    Flag_ecalBadCalibFilter
                    """
        rdf = (rdf.Filter(flag))
        return rdf
        
    def __loadRDF(self, filelist : List[str], treename : str = "Events", nFiles : int = -1, local : bool = False) -> RNode:
        if nFiles > 0:
            if len(filelist) < nFiles:
                raise ValueError("The filelist does not contain enough files")
            else:
                filelist = filelist[:nFiles]
        else:
            nFiles = len(filelist)
            
        if not local:
            filelist = ["root://cms-xrd-global.cern.ch/" + file for file in filelist[:nFiles]]
        else:
            filelist = filelist[:nFiles]
    
        self.chain = ROOT.TChain(treename)
        for file in filelist[:nFiles]:
            self.chain.Add(file)
        
        rdf = RDataFrame(self.chain)
        
        return rdf
 
    def __redo_JEC(self, jec : JEC_corrections) -> RNode:
        ROOT.gInterpreter.Declare('#include "src/JECRDF_code.h"')
        # Remember that this might alter the leading jets!
        ROOT.init_JEC(L1 = jec.L1, L2Relative = jec.L2Relative, L2L3 = jec.L2L3, nThreads = self.nThreads)

        rdf = (self.rdf.Define("JEC", "getJEC(rdfslot_, Jet_pt, Jet_eta, Jet_area, Rho_fixedGridRhoFastjetAll)")
               .Redefine("Jet_pt", "Jet_rawPt * JEC")
               .Redefine("Jet_rawFactor", "1.0 - 1.0 / JEC")
               .Redefine("Jet_order", "ROOT::VecOps::Argsort(Jet_pt)")
        )
        self.histograms["all"].extend([
            rdf.Histo1D(("JEC", "JEC;JEC;N_{events}", self.bins["response"]["n"], self.bins["response"]["pt"]), "JEC", "weight"),
            ])

        return rdf
    
    def __do_cut_golden_json(self, json_file : str) -> RNode:
        ROOT.gInterpreter.Declare('#include "src/JSONRDF_code.h"')
        ROOT.init_json(json_file)
        rdf = self.rdf.Filter("isGoodLumi(run, luminosityBlock)", "JSON Filter")
        return rdf
    
    def __do_cut_veto_map(self, veto_map_file : str) -> RNode:
        ROOT.gInterpreter.Declare('#include "src/JETVETOMAPS_code.h"')
        ROOT.init_vetomap(veto_map_file)
        rdf = self.rdf.Redefine("Jet_passesVetomap","isGoodVeto(Jet_eta, Jet_phi)")
        return rdf
    
    def do_smear_JER(self):
        # TODO: Implement this in C++
        pass
            
    def get_histograms(self) -> dict:
        print("Returning histograms")
        return self.histograms

    def run_histograms(self) -> "RDFAnalyzer":
        if not self.has_run:
            for trigger in self.trigger_list:
                print("Running histograms for trigger", trigger)
                RunGraphs(self.histograms[trigger])
            
            self.has_run = True
        else:
            print("Histograms have already been run")
        return self
            
    def do_inclusive(self) -> "RDFAnalyzer":
        # Create the inclusive histograms
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = (self.Flag_cut(rdf)).Redefine("Jet_pt", "Jet_pt[Jet_jetId >= 4]").Redefine("Jet_eta", "Jet_eta[Jet_jetId >= 4]")
            
            # Eta binned rdfs for pT distribution of jets
            eta_bins_for_pt = [(0.0, 1.3), (0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5),
                                 (2.5, 3.0), (3.0, 3.5), (3.5, 4.0), (4.0, 4.5), (4.5, 5.0)]
            eta_binned_rdfs = {}
            for i, val in enumerate(eta_bins_for_pt):
                eta_binned_rdfs[i] = (selected_rdf.Redefine("Jet_pt", f"Jet_pt[abs(Jet_eta) > {val[0]} && abs(Jet_eta) < {val[1]}]"))
                
            print("Creating inclusive histograms for trigger", trigger)
            self.histograms[trigger].extend([
                all_rdf.Histo2D(("Inclusive_EtaVsPt_all", "Inclusive_EtaVsPt;|#eta|;p_{T} (GeV);",
                                self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Jet_eta", "Jet_pt", "weight"),
                selected_rdf.Histo2D(("Inclusive_EtaVsPt_selected", "Inclusive_EtaVsPt;|#eta_{jet}|;p_{T} (GeV)",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                    "Jet_eta", "Jet_pt", "weight")
            ]
            )
            
            self.histograms[trigger].extend([
                eta_binned_rdfs[i].Histo1D(("Inclusive_Pt_eta_" + str(eta_bins_for_pt[i][1]), "Inclusive_pT_eta_" + str(eta_bins_for_pt[i][1]) + ";p_{T} (GeV)", 
                                            self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                           "Jet_pt", "weight") for i in range(len(eta_binned_rdfs))])
            
        return self
    
    def do_PFComposition(self) -> "RDFAnalyzer":
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = ((self.Flag_cut(rdf))
                            .Redefine("Jet_pt", "Jet_pt[Jet_jetId >= 4]")
                            .Redefine("Jet_eta", "Jet_eta[Jet_jetId >= 4]")
                            .Redefine("Jet_neHEF", "Jet_neHEF[Jet_jetId >= 4]")
                            .Redefine("Jet_neEmEF", "Jet_neEmEF[Jet_jetId >= 4]")
                            .Redefine("Jet_chHEF", "Jet_chHEF[Jet_jetId >= 4]")
                            .Redefine("Jet_chEmEF", "Jet_chEmEF[Jet_jetId >= 4]")
                            .Redefine("Jet_muEF", "Jet_muEF[Jet_jetId >= 4]")
                            )
            
            print("Creating PFComposition histograms for trigger", trigger)
            # TODO: This kind of behaviour of repeating histogram creation could be optimized
            self.histograms[trigger].extend([
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfilePt_all", "PFComposition_EtaVsPtVsProfilePt;|#eta|;p_{T} (GeV);p_{T} (GeV);", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_pt", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileRho_all", "PFComposition_EtaVsPtVsProfileRho;|#eta|;p_{T} (GeV);#rho;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Rho_fixedGridRhoFastjetAll", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNHF_all", "PFComposition_EtaVsPtVsProfileNHF;|#eta|;p_{T} (GeV);NHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_neHEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNEF_all", "PFComposition_EtaVsPtVsProfileNEF;|#eta|;p_{T} (GeV);NEF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_neEmEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCHF_all", "PFComposition_EtaVsPtVsProfileCHF;|#eta|;p_{T} (GeV);CHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_chHEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCEF_all", "PFComposition_EtaVsPtVsProfileCEF;|#eta|;p_{T} (GeV);CEF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_chEmEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileMUF_all", "PFComposition_EtaVsPtVsProfileMUF;|#eta|;p_{T} (GeV);MUF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_muEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfilePt_selected", "PFComposition_EtaVsPtVsProfilePt|#eta_{jet}|;p_{T} (GeV);p_{T} (GeV);", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_pt", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileRho_selected", "PFComposition_EtaVsPtVsProfileRho|#eta_{jet}|;p_{T} (GeV);#rho;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Rho_fixedGridRhoFastjetAll", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNHF_selected", "PFComposition_EtaVsPtVsProfileNHF|#eta_{jet}|;p_{T} (GeV);NHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_neHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNEF_selected", "PFComposition_EtaVsPtVsProfileNEF|#eta_{jet}|;p_{T} (GeV);NEF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_neEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCHF_selected", "PFComposition_EtaVsPtVsProfileCHF|#eta_{jet}|;p_{T} (GeV);CHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_chHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCEF_selected", "PFComposition_EtaVsPtVsProfileCEF|#eta_{jet}|;p_{T} (GeV);CEF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_chEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileMUF_selected", "PFComposition_EtaVsPtVsProfileMUF|#eta_{jet}|;p_{T} (GeV);MUF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_muEF", "weight")
            ])
            
        return self
    
    def __sample_cut(self, rdf : RNode):
        pass




