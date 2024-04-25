import ROOT
import correctionlib
from dataclasses import dataclass
from typing import List, Tuple, Dict
from RDFHelpers import get_bins, get_fill_range, find_era, update_run_bins
import numpy as np
import copy
correctionlib.register_pyroot_binding()


RDataFrame = ROOT.RDataFrame
RNode = ROOT.RDF.RNode
RunGraphs = ROOT.RDF.RunGraphs
initialize = ROOT.RDF.Experimental.Distributed.initialize


class RDFAnalyzer:
    def __init__(self, filelist : List[str],
                triggers : Dict = {},
                json_file : str = "",
                nFiles : int = -1,
                JEC : Dict = {},
                nThreads : int = 1,
                progress_bar : bool = False,
                isMC : bool = False,
                local : bool = False,
                system : str = "standard",
                run_raw : bool = False,
                selection_only : bool = True,
                header_dir : str = "src"
                ):

        self.nThreads = nThreads
        self.triggers = triggers
        self.trigger_list = copy.deepcopy(list(triggers.keys())) # Is deepcopy required if list(.keys())?
        self.histograms = {"all" : []} # format : {trigger : [histograms]}
        self.trigger_rdfs = {} # format : {trigger : rdf}. self.rdf is not initialized here due to order of operations
        self.has_run = False # possibly unnecessary
        self.chain = None
        self.isMC = isMC
        self.era = ""
        self.system = system
        self.run_raw= run_raw
        self.selection_only = selection_only
        self.header_dir = header_dir

        if not self.isMC:
            # Find the Run from filename
            self.era = find_era(filelist)
            frange = get_fill_range(self.era)
            self.bins = get_bins(fill_range=frange)
            print(f"Using fill range {frange} for run {self.era}")
        else:
            self.bins = get_bins()
        
        self.rdf = self.__loadRDF(filelist, nFiles = nFiles, local = local)
        if progress_bar:
            ROOT.RDF.Experimental.AddProgressBar(self.rdf)
            
        trigger_string = ""
        for trigger in triggers:
            trigger_string += triggers[trigger] + " || "
        trigger_string = trigger_string[:-4]

        if trigger_string == "":
            trigger_string = "1"
        
        self.trigger_list.append("all_trigs")

        # Initial variables
        self.rdf = (self.rdf.Define("weight", "genWeight" if self.isMC else "1.0")
                    .Filter("nJet > 0", "Only events with jets")
                    .Define("Jet_order", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Jet_pt))")
                    .Define("Jet_rawPt", "Jet_pt * (1.0-Jet_rawFactor)")
                    .Define("Jet_passesVetomap", "ROOT::VecOps::RVec<int>(Jet_pt.size(), 1)")
                    .Define("RawPuppiMET_polar", "ROOT::Math::Polar2DVectorF(RawPuppiMET_pt, RawPuppiMET_phi)")
                    .Define("all_trigs", trigger_string)
                    .Define("all", "1")
                    .Define("Jet_pt_leading", "Jet_pt[Jet_order[0]]")
                    .Define("Jet_eta_leading", "Jet_eta[Jet_order[0]]")
                    .Define("JEC", "1.0")
                    )
        if not self.isMC:
            self.bins = update_run_bins(self.rdf, self.bins)
        
        # MC cuts, to be implemented elsewhere
        if self.isMC:
            self.rdf = (self.rdf.Filter("fabs(PV_z - GenVtx_z) < 0.2", "Vertex_z_cut")
                        .Define("JetWithGen_eta", "Jet_eta[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_pt", "Jet_pt[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_phi", "Jet_phi[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_mass", "Jet_mass[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_rawFactor", "Jet_rawFactor[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_rawPt", "Jet_rawPt[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_area", "Jet_area[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_jetId", "Jet_jetId[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_order", "ROOT::VecOps::Argsort(JetWithGen_pt)")
                        .Define("JetWithGen_pt_leading", "Jet_pt[JetWithGen_order[0]]")
                        .Define("JetWithGen_neHEF", "Jet_neHEF[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_neEmEF", "Jet_neEmEF[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_chHEF", "Jet_chHEF[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_chEmEF", "Jet_chEmEF[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_muEF", "Jet_muEF[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_genJetIdx", "Jet_genJetIdx[Jet_genJetIdx >= 0]")
                        .Define("JetWithGen_passesVetomap", "ROOT::VecOps::RVec<int>(JetWithGen_pt.size(), 1)")
                        .Define("JetNoGen_pt", "Jet_pt[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_eta", "Jet_eta[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_phi", "Jet_phi[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_mass", "Jet_mass[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_rawFactor", "Jet_rawFactor[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_rawPt", "Jet_rawPt[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_area", "Jet_area[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_jetId", "Jet_jetId[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_order", "ROOT::VecOps::Argsort(JetNoGen_pt)")
                        .Define("JetNoGen_pt_leading", "Jet_pt[JetNoGen_order[0]]")
                        .Define("JetNoGen_neHEF", "Jet_neHEF[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_neEmEF", "Jet_neEmEF[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_chHEF", "Jet_chHEF[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_chEmEF", "Jet_chEmEF[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_muEF", "Jet_muEF[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_genJetIdx", "Jet_genJetIdx[Jet_genJetIdx < 0]")
                        .Define("JetNoGen_passesVetomap", "ROOT::VecOps::RVec<int>(JetNoGen_pt.size(), 1)")
                        .Filter("JetWithGen_pt.size() > 0", "At least one gen jet")
                        .Redefine("GenJet_pt", "ROOT::VecOps::Take(GenJet_pt, JetWithGen_genJetIdx)")
                        .Redefine("GenJet_eta", "ROOT::VecOps::Take(GenJet_eta, JetWithGen_genJetIdx)")
                        .Redefine("GenJet_phi", "ROOT::VecOps::Take(GenJet_phi, JetWithGen_genJetIdx)")
                        .Redefine("GenJet_mass", "ROOT::VecOps::Take(GenJet_mass, JetWithGen_genJetIdx)")
                        .Redefine("GenJet_partonFlavour", "ROOT::VecOps::Take(GenJet_partonFlavour, JetWithGen_genJetIdx)")
                        )
            
        # JECs and JERs
        if JEC:
            self.rdf = self.apply_correction(self.rdf, JEC)

        # Check if raw pt run
        if self.run_raw:
            self.rdf = (self.rdf.Redefine("Jet_pt", "Jet_pt * (1.0 - Jet_rawFactor)")
                        .Redefine("Jet_order", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Jet_pt))")
                        .Redefine("Jet_pt_leading", "Jet_pt[Jet_order[0]]")
            )

        # JSON cut
        if (json_file != "") and (not self.isMC):
            self.rdf = self.__do_cut_golden_json(json_file)

        # Triggers 
        for trigger in self.trigger_list:
            if trigger != "":
                self.trigger_rdfs[trigger] = self.rdf.Filter(trigger)
                self.histograms[trigger] = []

        # Only after JECs, common cuts etc. have been applied, we can create the inclusive rdf
        self.trigger_rdfs["all"] = self.rdf
        self.trigger_list.append("all")
        
        # Set gRandom seed
        ROOT.gRandom.SetSeed(12345)
        
    def do_inclusive_control(self) -> "RDFAnalyzer":
        print(f"Creating inclusive control histograms for {self.system}")
        filter_flags = (
            "Flag_goodVertices",
            "Flag_globalSuperTightHalo2016Filter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter",
            "Flag_BadPFMuonDzFilter",
            "Flag_hfNoisyHitsFilter",
            "Flag_eeBadScFilter",
            "Flag_ecalBadCalibFilter"
        )
        for flag in filter_flags:
            self.histograms["all"].extend([
                self.rdf.Histo1D((f"Control_{flag}", flag+";"+flag+";N_{events}", 2, 0, 2), flag, "weight")
            ])
        self.histograms["all"].extend([
            self.rdf.Histo1D((f"Control_jetID", "Jet_jetID;Jet_jetID;N_{events}", 8, 0, 8), "Jet_jetId", "weight"),
        ])

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
 
    def apply_correction(self, in_rdf: RNode, JEC: Dict) -> RNode:
        out_rdf = in_rdf 
        for cset in JEC.keys():
            file_name = JEC[cset]["file"]
            func_call = JEC[cset]["func_call"]
            rdf_call = JEC[cset]["rdf_call"]
            eval_call = JEC[cset]["eval_call"]
            
            ROOT.gInterpreter.Declare(f"""
            #ifndef INCLUDE_CORRECTIONLIB_CODE{cset}
            #define INCLUDE_CORRECTIONLIB_CODE{cset}
            auto correctionSet{cset} = correction::CorrectionSet::from_file("{file_name}");
            auto correctionEvaluator{cset} = correctionSet{cset}->at("{cset}");
            #endif
            """)
            ROOT.gInterpreter.Declare(f"""
            #ifndef INCLUDE_CORRECTIONLIB_FUNC{cset}
            #define INCLUDE_CORRECTIONLIB_FUNC{cset}
            """+func_call+""" {
            ROOT::VecOps::RVec<float> JECs;
                for (int i = 0; i < Jet_pt.size(); i++) {
                    JECs.push_back(correctionEvaluator"""+cset+"""->evaluate("""+eval_call+"""));
                }
                return JECs;
            }
            #endif
                                    """
            )

            out_rdf =( out_rdf.Redefine("JEC", "JEC * "+rdf_call)
                    .Redefine("Jet_pt", "Jet_rawPt * JEC")
                    .Redefine("Jet_rawFactor", "1.0 - 1.0 / JEC")
                    .Redefine("Jet_order", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Jet_pt))") 
                    .Redefine("Jet_pt_leading", "Jet_pt[Jet_order[0]]")
            )

        return out_rdf
    
    def __do_cut_golden_json(self, json_file : str) -> RNode:
        print("Applying golden JSON cut")
        print("JSON file:", json_file)

        ROOT.gInterpreter.Declare(f'''
        #ifndef INCLUDE_JSONRDF_CODE
        #define INCLUDE_JSONRDF_CODE
        #include "{self.header_dir}/JSONRDF_code.h"
        #endif
        ''')

        ROOT.init_json(json_file)
        rdf = self.rdf.Define("goldenJSON", "isGoodLumi(run, luminosityBlock)").Filter("goldenJSON", "JSON Filter")
        self.histograms["all"].extend([
            rdf.Histo1D(("GoldenJSON", "GoldenJSON;GoldenJSON;N_{events}", 2, 0, 2), "goldenJSON", "weight")
        ])
        return rdf
    
    def __do_cut_veto_map(self, veto_map_file : str) -> RNode:
        ROOT.gInterpreter.Declare(f'#include "{self.header_dir}/JETVETOMAPS_code.h"')
        ROOT.init_vetomap(veto_map_file)
        rdf = self.rdf.Redefine("Jet_passesVetomap","isGoodVeto(Jet_eta, Jet_phi)")

        if self.isMC:
            rdf = rdf.Redefine("JetWithGen_passesVetomap","isGoodVeto(JetWithGen_eta, JetWithGen_phi)")
            rdf = rdf.Redefine("JetNoGen_passesVetomap","isGoodVeto(JetNoGen_eta, JetNoGen_phi)")
            
        self.histograms["all"].extend([
            rdf.Histo1D(("VetoMap_VetoMap", "VetoMap;VetoMap;N_{events}", 2, 0, 2), "Jet_passesVetomap", "weight"),
        ])
        if self.isMC:
            self.histograms["all"].extend([
                rdf.Histo1D(("VetoMap_VetoMapWithGen", "VetoMapWithGen;VetoMapWithGen;N_{events}", 2, 0, 2), "JetWithGen_passesVetomap", "weight"),
                rdf.Histo1D(("VetoMap_VetoMapNoGen", "VetoMapNoGen;VetoMapNoGen;N_{events}", 2, 0, 2), "JetNoGen_passesVetomap", "weight")
            ])
        return rdf
    
    def get_histograms(self) -> dict:
        if not self.has_run:
            print("Please run histograms first")
            return {}
        else:
            return self.histograms

    def run_histograms(self) -> "RDFAnalyzer":
        if not self.has_run:
            print(f"Running histograms for system: {self.system}")
            for trigger in self.trigger_list:
                if trigger != "":
                    RunGraphs(self.histograms[trigger])
            
            self.has_run = True
        else:
            print("Histograms have already been run")
        return self
            

    def do_inclusive(self) -> "RDFAnalyzer":
        # Create the inclusive histograms
        # HOX! How do you know that the HLT jet isn't cut away with Jet_jetId or eta cut?
        print(f"Creating inclusive histograms for {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = ((self.Flag_cut(rdf))
                            .Redefine("Jet_pt_leading", "Jet_pt[Jet_jetId >= 4][0]")
                            .Redefine("Jet_eta_leading", "Jet_eta[Jet_jetId >= 4][0]")
                            .Redefine("Jet_pt", "Jet_pt[Jet_jetId >= 4]")
                            .Redefine("Jet_eta", "Jet_eta[Jet_jetId >= 4]")
                            .Redefine("Jet_phi", "Jet_phi[Jet_jetId >= 4]")
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
                
            if not self.selection_only:
                self.histograms[trigger].extend([
                    all_rdf.Histo2D(("Inclusive_EtaVsPt_all", "Inclusive_EtaVsPt;#eta;p_{T} (GeV);",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Jet_eta", "Jet_pt", "weight"),
                    all_rdf.Histo2D(("Inclusive_EtaVsPtlead_all", "Inclusive_EtaVsPtlead;#eta;p_{T,lead} (GeV);",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Jet_eta_leading", "Jet_pt_leading", "weight"),
                    all_rdf.Histo1D(("Inclusive_Pt_all", "Inclusive_pT_all;p_{T} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt", "weight"),
                    all_rdf.Histo1D(("Inclusive_Ptlead_all", "Inclusive_pTlead_all;p_{T,lead} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt_leading", "weight"),
                    all_rdf.Histo2D(("Inclusive_EtaVsPhi_all", "Inclusive_EtaVsPhi;#eta;#phi;",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                    "Jet_eta", "Jet_phi", "weight"),
                ])
            self.histograms[trigger].extend([
                selected_rdf.Histo2D(("Inclusive_EtaVsPt_selected", "Inclusive_EtaVsPt;#eta_{jet};p_{T} (GeV)",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                    "Jet_eta", "Jet_pt", "weight"),
                selected_rdf.Histo2D(("Inclusive_EtaVsPtlead_selected", "Inclusive_EtaVsPtlead;#eta_{jet};p_{T,lead} (GeV)",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                                    "Jet_eta_leading", "Jet_pt_leading", "weight"),
                selected_rdf.Histo1D(("Inclusive_Pt_selected", "Inclusive_pT_selected;p_{T} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt", "weight"),
                selected_rdf.Histo1D(("Inclusive_Ptlead_selected", "Inclusive_pTlead_selected;p_{T,lead} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]), "Jet_pt_leading", "weight"),
                selected_rdf.Histo2D(("Inclusive_EtaVsPhi_selected", "Inclusive_EtaVsPhi;#eta;#phi;",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                    "Jet_eta", "Jet_phi", "weight")
            
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
    
    def do_MC(self) -> "RDFAnalyzer":
        if not self.isMC:
            raise ValueError("Trying to use MC specific functions on data")
        
        print(f"Creating MC histograms for system: {self.system}")
        
        for trigger, rdf in self.trigger_rdfs.items():
            rdf = self._MC_cut(rdf)
            all_rdf = rdf
            selected_rdf = ((self.Flag_cut(rdf)))
            
            for inner_rdf, rdf_name in zip([all_rdf, selected_rdf], ["all", "selected"]):
                
                self.histograms[trigger].extend([
                    inner_rdf.Histo1D((f"MC_GenJet_pt_{rdf_name}", "GenJet_pt;p_{T} (GeV);N_{events}", self.bins["pt"]["n"], self.bins["pt"]["bins"]), "GenJet_pt", "weight"),
                    inner_rdf.Histo1D((f"MC_GenJet_eta_{rdf_name}", "GenJet_eta;|#eta|;N_{events}", self.bins["eta"]["n"], self.bins["eta"]["bins"]), "GenJet_eta", "weight"),
                    inner_rdf.Histo1D((f"MC_GenJet_phi_{rdf_name}", "GenJet_phi;#phi;N_{events}", self.bins["phi"]["n"], self.bins["phi"]["bins"]), "GenJet_phi", "weight"),
                    inner_rdf.Histo1D((f"MC_GenJet_partonFlavour_{rdf_name}", "GenJet_partonFlavour;Parton flavour;N_{events}", 44, -22, 22), "GenJet_partonFlavour", "weight"),
                    inner_rdf.Histo1D((f"MC_Response_{rdf_name}", "Response;Response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), "Response", "weight"),
                    inner_rdf.Histo1D((f"MC_Response_leading_{rdf_name}", "Response_leading;Response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_leading", "weight"),
                    inner_rdf.Histo1D((f"MC_Response_twoLeading_{rdf_name}", "Response_twoLeading;Response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_twoLeading", "weight"),
                    inner_rdf.Histo1D((f"MC_Response_threeLeading_{rdf_name}", "Response_threeLeading;Response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_threeLeading", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsPt_{rdf_name}", "PtVsResponse;GenJet p_{T} (GeV);Response;", self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "GenJet_pt", "Response", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsEta_{rdf_name}", "EtaVsResponse;GenJet |#eta|;Response;", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "GenJet_eta", "Response", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsPt_lead_{rdf_name}", "PtVsResponse;GenJet p_{T} (GeV);Response;", self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_leadingPt", "Response_leading", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsEta_lead_{rdf_name}", "EtaVsResponse;GenJet |#eta|;Response;", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_leadingEta", "Response_leading", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsPt_twoLead_{rdf_name}", "PtVsResponse;GenJet p_{T} (GeV);Response;", self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_twoLeadingPts", "Response_twoLeading", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsEta_twoLead_{rdf_name}", "EtaVsResponse;GenJet |#eta|;Response;", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_twoLeadingEtas", "Response_twoLeading", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsPt_threeLead_{rdf_name}", "PtVsResponse;GenJet p_{T} (GeV);Response;", self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_threeLeadingPts", "Response_threeLeading", "weight"),
                    inner_rdf.Histo2D((f"MC_ResponseVsEta_threeLead_{rdf_name}", "EtaVsResponse;GenJet |#eta|;Response;", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), "Response_threeLeadingEtas", "Response_threeLeading", "weight"),
                ])
                
        return self
                
    def _MC_cut(self, rdf : RNode) -> RNode:
        # TODO: Account for changes in order of jets
        rdf_MC = (rdf.Filter("nGenJet > 0", "GenJet_cut")
                    .Define("Response", "JetWithGen_pt / GenJet_pt")
                    .Define("Response_leading", "JetWithGen_pt[0] / GenJet_pt[0]")
                    .Define("Response_twoLeading", "JetWithGen_pt.size() > 1 ? ROOT::VecOps::Take(JetWithGen_pt, 2) / ROOT::VecOps::Take(GenJet_pt, 2) : ROOT::VecOps::Take(JetWithGen_pt, 1) / ROOT::VecOps::Take(GenJet_pt, 1)")
                    .Define("Response_threeLeading", "JetWithGen_pt.size() > 2 ? ROOT::VecOps::Take(JetWithGen_pt, 3) / ROOT::VecOps::Take(GenJet_pt, 3) : Response_twoLeading")
                    .Define("Response_leadingPt", "ROOT::VecOps::Take(GenJet_pt, 1)")
                    .Define("Response_twoLeadingPts", "JetWithGen_pt.size() > 1 ? ROOT::VecOps::Take(GenJet_pt, 2) : ROOT::VecOps::Take(GenJet_pt, 1)")
                    .Define("Response_threeLeadingPts", "JetWithGen_pt.size() > 2 ? ROOT::VecOps::Take(GenJet_pt, 3) : Response_twoLeadingPts")
                    .Define("Response_leadingEta", "JetWithGen_eta[0]")
                    .Define("Response_twoLeadingEtas", "JetWithGen_pt.size() > 1 ? ROOT::VecOps::Take(GenJet_eta, 2) : ROOT::VecOps::Take(GenJet_eta, 1)")
                    .Define("Response_threeLeadingEtas", "JetWithGen_pt.size() > 2 ? ROOT::VecOps::Take(GenJet_eta, 3) : Response_twoLeadingEtas")
                    .Define("Leading_uds", "abs(GenJet_partonFlavour[0]) == 1 || abs(GenJet_partonFlavour[0]) == 2 || abs(GenJet_partonFlavour[0]) == 3")
                    .Define("Leading_ud", "abs(GenJet_partonFlavour[0]) == 1 || abs(GenJet_partonFlavour[0]) == 2")
                    .Define("Leading_s", "abs(GenJet_partonFlavour[0]) == 3")
                    .Define("Leading_c", "abs(GenJet_partonFlavour[0]) == 4")
                    .Define("Leading_b", "abs(GenJet_partonFlavour[0]) == 5")
                    .Define("Leading_g", "abs(GenJet_partonFlavour[0]) == 21")
        )
        return rdf_MC
    
    def do_RunsAndLumis(self) -> "RDFAnalyzer":
        print(f"Creating run and lumi histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = (self.Flag_cut(rdf))
            if not self.selection_only:
                self.histograms[trigger].extend([
                    all_rdf.Histo1D(("RunAndLumi_Run_all", "Run_all;Run;N_{events}", self.bins["runs"]["n"], self.bins["runs"]["bins"]), "run", "weight"),
                    all_rdf.Histo1D(("RunAndLumi_Lumi_all", "Lumi_all;Lumi;N_{events}", self.bins["lumi"]["n"], self.bins["lumi"]["bins"]), "luminosityBlock", "weight"),
                    all_rdf.Histo1D(("RunAndLumi_BunchCrossing_all", "BunchCrossing_all;BunchCrossing;N_{events}", self.bins["bx"]["n"], self.bins["bx"]["bins"]), "bunchCrossing", "weight"),
                    all_rdf.Histo2D(("RunAndLumi_RunVsBunchCrossing_all", "RunVsBunchCrossing_all;Run;BunchCrossing;N_{events}", self.bins["runs"]["n"], self.bins["runs"]["bins"], self.bins["bx"]["n"], self.bins["bx"]["bins"]), "run", "bunchCrossing", "weight"),
                ])
            self.histograms[trigger].extend([
                selected_rdf.Histo1D(("RunAndLumi_Run_selected", "Run_selected;Run;N_{events}", self.bins["runs"]["n"], self.bins["runs"]["bins"]), "run", "weight"),
                selected_rdf.Histo1D(("RunAndLumi_Lumi_selected", "Lumi_selected;Lumi;N_{events}", self.bins["lumi"]["n"], self.bins["lumi"]["bins"]), "luminosityBlock", "weight"),
                selected_rdf.Histo1D(("RunAndLumi_BunchCrossing_selected", "BunchCrossing_selected;BunchCrossing;N_{events}", self.bins["bx"]["n"], self.bins["bx"]["bins"]), "bunchCrossing", "weight"),
                selected_rdf.Histo2D(("RunAndLumi_RunVsBunchCrossing_selected", "RunVsBunchCrossing_selected;Run;BunchCrossing;N_{events}", self.bins["runs"]["n"], self.bins["runs"]["bins"], self.bins["bx"]["n"], self.bins["bx"]["bins"]), "run", "bunchCrossing", "weight")
            ])
        return self
    
    def do_PFComposition(self) -> "RDFAnalyzer":
        print("Creating PFComposition histograms for system: ", self.system)
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
            if not self.selection_only:
                self.histograms[trigger].extend([
                    all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileRho_all", "PFComposition_EtaVsPtVsProfileRho;#eta;p_{T} (GeV);#rho;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                        "Jet_eta", "Jet_pt", "Rho_fixedGridRhoFastjetAll", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNHF_all", "PFComposition_EtaVsPtVsProfileNHF;#eta;p_{T} (GeV);NHF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                        "Jet_eta", "Jet_pt", "Jet_neHEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNEF_all", "PFComposition_EtaVsPtVsProfileNEF;#eta;p_{T} (GeV);NEF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                        "Jet_eta", "Jet_pt", "Jet_neEmEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCHF_all", "PFComposition_EtaVsPtVsProfileCHF;#eta;p_{T} (GeV);CHF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                        "Jet_eta", "Jet_pt", "Jet_chHEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCEF_all", "PFComposition_EtaVsPtVsProfileCEF;#eta;p_{T} (GeV);CEF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                        "Jet_eta", "Jet_pt", "Jet_chEmEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileMUF_all", "PFComposition_EtaVsPtVsProfileMUF;#eta;p_{T} (GeV);MUF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                        "Jet_eta", "Jet_pt", "Jet_muEF", "weight"),
                ])
            # TODO: This kind of behaviour of repeating histogram creation could be optimized
            self.histograms[trigger].extend([
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileRho_selected", "PFComposition_EtaVsPtVsProfileRho;#eta_{jet};p_{T} (GeV);#rho;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Rho_fixedGridRhoFastjetAll", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNHF_selected", "PFComposition_EtaVsPtVsProfileNHF;#eta_{jet};p_{T} (GeV);NHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_neHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNEF_selected", "PFComposition_EtaVsPtVsProfileNEF;#eta_{jet};p_{T} (GeV);NEF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_neEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCHF_selected", "PFComposition_EtaVsPtVsProfileCHF;#eta_{jet};p_{T} (GeV);CHF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_chHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCEF_selected", "PFComposition_EtaVsPtVsProfileCEF;#eta_{jet};p_{T} (GeV);CEF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_chEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileMUF_selected", "PFComposition_EtaVsPtVsProfileMUF;#eta_{jet};p_{T} (GeV);MUF;", 
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]), 
                    "Jet_eta", "Jet_pt", "Jet_muEF", "weight")
            ])
            
            if not self.selection_only:
                self.histograms[trigger].extend([
                    all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfilePt_all", "PFComposition_EtaVsPhiVsProfilePt;#eta;#phi;p_{T} (GeV);", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                        "Jet_eta", "Jet_phi", "Jet_pt", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileRho_all", "PFComposition_EtaVsPhiVsProfileRho;#eta;#phi;#rho;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                        "Jet_eta", "Jet_phi", "Rho_fixedGridRhoFastjetAll", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNHF_all", "PFComposition_EtaVsPhiVsProfileNHF;#eta;#phi;NHF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                        "Jet_eta", "Jet_phi", "Jet_neHEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNEF_all", "PFComposition_EtaVsPhiVsProfileNEF;#eta;#phi;NEF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                        "Jet_eta", "Jet_phi", "Jet_neEmEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCHF_all", "PFComposition_EtaVsPhiVsProfileCHF;#eta;#phi;CHF;", 
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                        "Jet_eta", "Jet_phi", "Jet_chHEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCEF_all", "PFComposition_EtaVsPhiVsProfileCEF;#eta;#phi;CEF;",
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                        "Jet_eta", "Jet_phi", "Jet_chEmEF", "weight"),
                    all_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileMUF_all", "PFComposition_EtaVsPhiVsProfileMUF;#eta;#phi;MUF;",
                        self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                        "Jet_eta", "Jet_phi", "Jet_muEF", "weight"),
                ])
            # Same plots as above but eta vs phi
            self.histograms[trigger].extend([
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfilePt_selected", "PFComposition_EtaVsPhiVsProfilePt;#eta_{jet};#phi;p_{T} (GeV);",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_pt", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileRho_selected", "PFComposition_EtaVsPhiVsProfileRho;#eta_{jet};#phi;#rho;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Rho_fixedGridRhoFastjetAll", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNHF_selected", "PFComposition_EtaVsPhiVsProfileNHF;#eta_{jet};#phi;NHF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_neHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileNEF_selected", "PFComposition_EtaVsPhiVsProfileNEF;#eta_{jet};#phi;NEF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_neEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCHF_selected", "PFComposition_EtaVsPhiVsProfileCHF;#eta_{jet};#phi;CHF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_chHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileCEF_selected", "PFComposition_EtaVsPhiVsProfileCEF;#eta_{jet};#phi;CEF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_chEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPhiVsProfileMUF_selected", "PFComposition_EtaVsPhiVsProfileMUF;#eta_{jet};#phi;MUF;",
                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                    "Jet_eta", "Jet_phi", "Jet_muEF", "weight")
            ])
            
        return self
    
    def __sample_cut(self, rdf : RNode):
        pass
    
    def do_sample_control(self, rdf : RNode):
        pass




