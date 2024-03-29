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
                local : bool = False,
                ) -> "RDFAnalyzer":
        self.nThreads = nThreads
        self.trigger_list = trigger_list
        self.histograms = {"all" : []} # format : {trigger : [histograms]}
        self.trigger_rdfs = {} # format : {trigger : rdf}. self.rdf is not initialized here due to order of operations
        self.has_run = False # possibly unnecessary
        self.chain = None
        self.bins = get_bins()
        self.JEC_included = False
        self.isMC = isMC
        print(isMC)
        
        self.rdf = self.__loadRDF(filelist, nFiles = nFiles, local = local)
        if progress_bar:
            ROOT.RDF.Experimental.AddProgressBar(self.rdf)
            
        trigger_string = " || ".join(trigger_list)
        # print(trigger_string)
        self.trigger_list.append("all_trigs")

        # Initial variables
        self.rdf = (self.rdf.Define("weight", "genWeight" if self.isMC else "1.0")
                    .Define("Jet_order", "ROOT::VecOps::Argsort(Jet_pt)")
                    .Define("Jet_rawPt", "Jet_pt * (1.0-Jet_rawFactor)")
                    .Define("Jet_passesVetomap", "ROOT::VecOps::RVec<int>(Jet_pt.size(), 1)")
                    .Define("RawPuppiMET_polar", "ROOT::Math::Polar2DVectorF(RawPuppiMET_pt, RawPuppiMET_phi)")
                    .Define("all_trigs", trigger_string)
                    .Define("Jet_pt_leading", "Jet_pt[Jet_order[0]]")
                    )
        
        # MC cuts, to be implemented elsewhere
        if self.isMC:
            self.rdf = (self.rdf.Filter("fabs(PV_z - GenVtx_z) < 0.2", "Vertex_z_cut")
                        .Redefine("Jet_pt", "Jet_pt[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_eta", "Jet_eta[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_phi", "Jet_phi[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_mass", "Jet_mass[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_rawFactor", "Jet_rawFactor[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_rawPt", "Jet_rawPt[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_area", "Jet_area[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_jetId", "Jet_jetId[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_order", "ROOT::VecOps::Argsort(Jet_pt)")
                        .Redefine("Jet_pt_leading", "Jet_pt[Jet_order[0]]")
                        .Redefine("Jet_neHEF", "Jet_neHEF[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_neEmEF", "Jet_neEmEF[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_chHEF", "Jet_chHEF[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_chEmEF", "Jet_chEmEF[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_muEF", "Jet_muEF[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_genJetIdx", "Jet_genJetIdx[Jet_genJetIdx >= 0]")
                        .Redefine("Jet_passesVetomap", "ROOT::VecOps::RVec<int>(Jet_pt.size(), 1)")
                        .Redefine("GenJet_pt", "ROOT::VecOps::Take(GenJet_pt, Jet_genJetIdx)")
                        .Redefine("GenJet_eta", "ROOT::VecOps::Take(GenJet_eta, Jet_genJetIdx)")
                        .Redefine("GenJet_phi", "ROOT::VecOps::Take(GenJet_phi, Jet_genJetIdx)")
                        .Redefine("GenJet_mass", "ROOT::VecOps::Take(GenJet_mass, Jet_genJetIdx)")
                        .Redefine("GenJet_partonFlavour", "ROOT::VecOps::Take(GenJet_partonFlavour, Jet_genJetIdx)")
                        )
            
        
        if not JEC.check_empty():
            print(JEC.L1, JEC.L2Relative, JEC.L2L3)
            # clean_JEC()
            # compile_JEC()
            load_JEC()
            self.rdf = self.__redo_JEC(JEC)

        if (json_file != "") and (not self.isMC):
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
        return rdf.Filter(flag, "Flag Filter")
        
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
        if not JEC_included:
            ROOT.gInterpreter.Declare('#include "src/JECRDF_code.h"')
        # Remember that this might alter the leading jets!
        ROOT.init_JEC(L1 = jec.L1, L2Relative = jec.L2Relative, L2L3 = jec.L2L3, nThreads = self.nThreads)

        rdf = (self.rdf.Define("JEC", "getJEC(rdfslot_, Jet_pt, Jet_eta, Jet_area, Rho_fixedGridRhoFastjetAll)")
               .Redefine("Jet_pt", "Jet_rawPt * JEC")
               .Redefine("Jet_rawFactor", "1.0 - 1.0 / JEC")
               .Redefine("Jet_order", "ROOT::VecOps::Argsort(Jet_pt)")
               .Redefine("Jet_pt_leading", "Jet_pt[Jet_order[0]]")
        )
        self.histograms["all"].extend([
            rdf.Histo1D(("JEC", "JEC;JEC;N_{events}", self.bins["response"]["n"], self.bins["response"]["pt"]), "JEC", "weight"),
            ])

        return rdf
    
    def __do_cut_golden_json(self, json_file : str) -> RNode:
        print("Applying golden JSON cut")
        print("JSON file:", json_file)
        ROOT.gInterpreter.Declare('#include "src/JSONRDF_code.h"')
        ROOT.init_json(json_file)
        rdf = self.rdf.Define("goldenJSON", "isGoodLumi(run, luminosityBlock)").Filter("goldenJSON", "JSON Filter")
        self.histograms["all"].extend([
            rdf.Histo1D(("GoldenJSON", "GoldenJSON;GoldenJSON;N_{events}", 2, 0, 2), "goldenJSON", "weight")
        ])
        return rdf
    
    def __do_cut_veto_map(self, veto_map_file : str) -> RNode:
        ROOT.gInterpreter.Declare('#include "src/JETVETOMAPS_code.h"')
        ROOT.init_vetomap(veto_map_file)
        rdf = self.rdf.Redefine("Jet_passesVetomap","isGoodVeto(Jet_eta, Jet_phi)")
        self.histograms["all"].extend([
            rdf.Histo1D(("VetoMap", "VetoMap;VetoMap;N_{events}", 2, 0, 2), "Jet_passesVetomap", "weight")
        ])
        return rdf
    
    def do_smear_JER(self):
        # TODO: Implement this in C++
        if not JEC_included:
            ROOT.gInterpreter.Declare('#include "src/JER_code.h"')
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
        # HOX! How do you know that the HLT jet isn't cut away with Jet_jetId or eta cut?
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = ((self.Flag_cut(rdf))
                            .Redefine("Jet_pt_leading", "Jet_pt[Jet_jetId >= 4][0]")
                            .Redefine("Jet_pt", "Jet_pt[Jet_jetId >= 4]")
                            .Redefine("Jet_eta", "Jet_eta[Jet_jetId >= 4]")
                            )
            
            # Eta binned rdfs for pT distribution of jets
            eta_bins_for_pt = [(0.0, 1.3), (0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5),
                                 (2.5, 3.0), (3.0, 3.5), (3.5, 4.0), (4.0, 4.5), (4.5, 5.0)]
            eta_binned_rdfs = {}
            eta_binned_lead_rdfs = {}
            for i, val in enumerate(eta_bins_for_pt):
                eta_binned_rdfs[i] = (selected_rdf.Redefine("Jet_pt", f"Jet_pt[abs(Jet_eta) > {val[0]} && abs(Jet_eta) < {val[1]}]")
                )
                eta_binned_lead_rdfs[i] = (selected_rdf.Filter(f"abs(Jet_eta[0]) > {val[0]} && abs(Jet_eta[0]) < {val[1]}")
                )
                
            print("Creating inclusive histograms for trigger", trigger)
            self.histograms[trigger].extend([
                all_rdf.Histo2D(("Inclusive_EtaVsPt_all", "Inclusive_EtaVsPt;|#eta|;p_{T} (GeV);",
                                self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Jet_eta", "Jet_pt", "weight"),
                all_rdf.Histo2D(("Inclusive_EtaVsPtlead_all", "Inclusive_EtaVsPtlead;|#eta|;p_{T,lead} (GeV);",
                                self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Jet_eta", "Jet_pt_leading", "weight"),
                all_rdf.Histo1D(("Inclusive_Pt_all", "Inclusive_pT_all;p_{T} (GeV);N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt", "weight"),
                all_rdf.Histo1D(("Inclusive_Ptlead_all", "Inclusive_pTlead_all;p_{T,lead} (GeV);N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt_leading", "weight"),
                selected_rdf.Histo2D(("Inclusive_EtaVsPt_selected", "Inclusive_EtaVsPt;|#eta_{jet}|;p_{T} (GeV)",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                    "Jet_eta", "Jet_pt", "weight"),
                selected_rdf.Histo2D(("Inclusive_EtaVsPtlead_selected", "Inclusive_EtaVsPtlead;|#eta_{jet}|;p_{T,lead} (GeV)",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                    "Jet_eta", "Jet_pt_leading", "weight"),
                selected_rdf.Histo1D(("Inclusive_Pt_selected", "Inclusive_pT_selected;p_{T} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt", "weight"),
                selected_rdf.Histo1D(("Inclusive_Ptlead_selected", "Inclusive_pTlead_selected;p_{T,lead} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt_leading", "weight")
            
            ]
            )
            
            self.histograms[trigger].extend([
                eta_binned_rdfs[i].Histo1D(("Inclusive_Pt_eta_" + str(eta_bins_for_pt[i][1]), "Inclusive_pT_eta_" + str(eta_bins_for_pt[i][1]) + ";p_{T} (GeV)", 
                                            self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                           "Jet_pt", "weight") for i in range(len(eta_binned_rdfs))])
            
            self.histograms[trigger].extend([
                eta_binned_lead_rdfs[i].Histo1D(("Inclusive_Ptlead_eta_" + str(eta_bins_for_pt[i][1]), "Inclusive_pTlead_eta_" + str(eta_bins_for_pt[i][1]) + ";p_{T,lead} (GeV)",
                                            self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                            "Jet_pt_leading", "weight") for i in range(len(eta_binned_lead_rdfs))])
            
        return self
    
    def do_RunsAndLumis(self) -> "RDFAnalyzer":
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = (self.Flag_cut(rdf))
            print("Creating run and lumi histograms for trigger", trigger)
            self.histograms[trigger].extend([
                all_rdf.Histo1D(("Run_all", "Run_all;Run;N_{events}", 373075-355100, 355100, 373075), "run", "weight"),
                all_rdf.Histo1D(("Lumi_all", "Lumi_all;Lumi;N_{events}", 1501, 0, 1500), "luminosityBlock", "weight"),
                selected_rdf.Histo1D(("Run_selected", "Run_selected;Run;N_{events}", 373075-355100, 355100, 373075), "run", "weight"),
                selected_rdf.Histo1D(("Lumi_selected", "Lumi_selected;Lumi;N_{events}", 1501, 0, 1500), "luminosityBlock", "weight")
            ])
        return self
    
    def do_PFComposition(self) -> "RDFAnalyzer":
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = ((self.Flag_cut(rdf))
                            .Redefine("Jet_pt", "Jet_pt[Jet_jetId >= 4]")
                            .Redefine("Jet_eta", "Jet_eta[Jet_jetId >= 4]")
                            .Redefine("Jet_phi", "Jet_phi[Jet_jetId >= 4]")
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
            
            # Same plots as above but eta vs phi
            self.histograms[trigger].extend([
                all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfilePt_all", "PFComposition_EtaVsPhiVsProfilePt;|#eta|;#phi;p_{T} (GeV);", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_pt", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileRho_all", "PFComposition_EtaVsPhiVsProfileRho;|#eta|;#phi;#rho;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Rho_fixedGridRhoFastjetAll", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNHF_all", "PFComposition_EtaVsPhiVsProfileNHF;|#eta|;#phi;NHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_neHEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNEF_all", "PFComposition_EtaVsPhiVsProfileNEF;|#eta|;#phi;NEF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_neEmEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCHF_all", "PFComposition_EtaVsPhiVsProfileCHF;|#eta|;#phi;CHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_chHEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCEF_all", "PFComposition_EtaVsPhiVsProfileCEF;|#eta|;#phi;CEF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_chEmEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileMUF_all", "PFComposition_EtaVsPhiVsProfileMUF;|#eta|;#phi;MUF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_muEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfilePt_selected", "PFComposition_EtaVsPhiVsProfilePt|#eta_{jet}|;#phi;p_{T} (GeV);",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_pt", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileRho_selected", "PFComposition_EtaVsPhiVsProfileRho|#eta_{jet}|;#phi;#rho;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Rho_fixedGridRhoFastjetAll", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNHF_selected", "PFComposition_EtaVsPhiVsProfileNHF|#eta_{jet}|;#phi;NHF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_neHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNEF_selected", "PFComposition_EtaVsPhiVsProfileNEF|#eta_{jet}|;#phi;NEF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_neEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCHF_selected", "PFComposition_EtaVsPhiVsProfileCHF|#eta_{jet}|;#phi;CHF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_chHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCEF_selected", "PFComposition_EtaVsPhiVsProfileCEF|#eta_{jet}|;#phi;CEF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_chEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileMUF_selected", "PFComposition_EtaVsPhiVsProfileMUF|#eta_{jet}|;#phi;MUF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_muEF", "weight")
            ])
            
        return self
    
    def __sample_cut(self, rdf : RNode):
        pass




