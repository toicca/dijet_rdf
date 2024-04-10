import ROOT
from typing import List
import numpy as np
from RDFAnalyzer import RDFAnalyzer, JEC_corrections, RNode
    
RDataFrame = ROOT.RDataFrame
RunGraphs = ROOT.RDF.RunGraphs
RNode = ROOT.RDF.RNode
    
class DijetAnalyzer(RDFAnalyzer):
    def __init__(self, filelist : List[str],
                trigger_list : List[str],
                json_file : str,
                nFiles : int = -1,
                JEC : JEC_corrections = JEC_corrections("", "", ""),
                nThreads : int = 1,
                progress_bar : bool = False,
                isMC : bool = False,
                local : bool = False
                ):
        super().__init__(filelist, trigger_list, json_file, nFiles, JEC, nThreads, progress_bar, isMC=isMC, local=local, system="dijet")
        
    def Flag_cut(self, rdf: RNode) -> RNode:
        return super().Flag_cut(rdf)
 

    def do_DB(self) -> "DijetAnalyzer":
        system = "dijet"
        print(f"Creating DB histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))
            pT_binLabels = ["average_Pt_dijet", "Jet_pt_tag", "Jet_pt_probe"]

            db_rdf = (db_rdf.Define(f"response_DB_{system}", f"(1.0 + Asymmetry_{system}) / (1.0 - Asymmetry_{system})")
                    )
            
            
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"DB_{system}_Response", "DB_"+ str(system) + "_response;response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), f"response_DB_{system}", "weight"),
                db_rdf.Histo2D((f"DB_{system}_EtaprobeVsResponse", "DB_"+ str(system) + "_EtaVsResponse;#eta_{probe};response",
                                self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                               f"Jet_eta_probe_{system}", f"response_DB_{system}", "weight"),
                db_rdf.Histo1D((f"DB_{system}_Asymmetry", "DB_"+ str(system) + "_Asymmetry;Asymmetry;N_{events}", self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]), 
                               f"Asymmetry_{system}", "weight"),
                db_rdf.Profile2D((f"DB_{system}_EtaprobeVsPhiprobeVsAsymmetry", "DB_"+ str(system) + "_EtaVsPhiVsAsymmetry;#eta_{probe};#phi_{probe};Asymmetry", 
                                self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                                f"Jet_eta_probe_{system}", f"Jet_phi_probe_{system}", f"Asymmetry_{system}", "weight"),
            ])
            for label in pT_binLabels:
                self.histograms[trigger].extend([
                    db_rdf.Histo2D((f"DB_{system}_PtVsResponse_{label}", "DB_"+str(system)+"_PtVsResponse_"+label+";p_{T, "+label+"} (GeV);response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["response"]["n"], self.bins["response"]["bins"]),
                                    label, f"response_DB_{system}", "weight"),
                    db_rdf.Histo3D((f"DB_{system}_PtVsEtaVsResponse_PtBin{label}", "MPF_"+str(system)+"_PtVsEtaVsResponse_"+label+";p_{T, "+label+"} (GeV);#eta_{probe};response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"],
                                    self.bins["response"]["n"], self.bins["response"]["bins"]),
                                    label, f"Jet_eta_probe_{system}", f"response_DB_{system}", "weight"),
                    db_rdf.Histo3D((f"DB_{system}_PtVsEtaVsAsymmetry_PtBin{label}", "MPF_"+str(system)+"_PtVsEtaVsAsymmetry_"+label+";p_{T, "+label+"} (GeV);#eta_{probe};response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"],
                                    self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                    label, f"Jet_eta_probe_{system}", f"Asymmetry_{system}", "weight"),
                    ]
                )

        return self
    
    def do_MPF(self) -> "DijetAnalyzer":
        system = "dijet"
        print(f"Creating MPF histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))
            pT_binLabels = ["average_Pt_dijet", "Jet_pt_tag", "Jet_pt_probe"]

            db_rdf = (db_rdf.Define(f"response_MPF_{system}", f"(1.0 + B_{system}) / (1.0 - B_{system})")
                    )
            
            
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"MPF_{system}_Response", "MPF_"+ str(system) + "_response;response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), 
                               f"response_MPF_{system}", "weight"),

                db_rdf.Histo2D((f"MPF_{system}_EtaprobeVsResponse", "DB_"+ str(system) + "_EtaprobeVsResponse;#eta_{probe};response", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), 
                               f"Jet_eta_probe_{system}", f"response_MPF_{system}", "weight"),
                db_rdf.Histo1D((f"MPF_{system}_B", "MPF_"+ str(system) + "_B;B;N_{events}", self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                 f"B_{system}", "weight"),
                db_rdf.Histo1D((f"MPF_{system}_avpResponse", "MPF_"+ str(system) + "_avpResponse;response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]),
                               "B_response", "weight"),
                db_rdf.Profile1D((f"MPF_{system}_PtVsAvpResponse", "MPF_"+ str(system) + "_PtVsAvpResponse;p_{T};response;N_{events}", self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "average_Pt_dijet", "B_response", "weight"),
                db_rdf.Profile1D((f"MPF_{system}_TagPtVsAvpResponse", "MPF_"+ str(system) + "_PtVsAvpResponse;p_{T};response;N_{events}", self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Jet_pt_tag", "B_response", "weight"),
                db_rdf.Profile1D((f"MPF_{system}_ProbePtVsAvpResponse", "MPF_"+ str(system) + "_PtVsAvpResponse;p_{T};response;N_{events}", self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Jet_pt_probe", "B_response", "weight"),
                db_rdf.Profile2D((f"MPF_{system}_EtaprobeVsPhiprobeVsB", "MPF_"+ str(system) + "_EtaprobeVsPhiprobeVsB;#eta_{probe};#phi_{probe};B", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                    f"Jet_eta_probe_{system}", f"Jet_phi_probe_{system}", f"B_{system}", "weight"),
                
            ])
            for label in pT_binLabels:
                self.histograms[trigger].extend([
                    db_rdf.Histo2D((f"MPF_{system}_PtVsResponse_{label}", "MPF_"+str(system)+"_PtVsResponse_"+label+";p_{T, "+label+"} (GeV);response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["response"]["n"], self.bins["response"]["bins"]),
                                    label, f"response_MPF_{system}", "weight"),
                    db_rdf.Histo3D((f"MPF_{system}_PtVsEtaVsResponse_PtBin{label}", "MPF_"+str(system)+"_PtVsEtaVsResponse_"+label+";p_{T, "+label+"} (GeV)|#eta_{probe};response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"],
                                    self.bins["response"]["n"], self.bins["response"]["bins"]),
                                    label, f"Jet_eta_probe_{system}", f"response_MPF_{system}", "weight"),
                    db_rdf.Histo3D((f"MPF_{system}_PtVsEtaVsB_PtBin{label}", "MPF_"+str(system)+"_PtVsEtaVsB_"+label+";p_{T, "+label+"} (GeV);#eta_{probe};response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"],
                                    self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                    label, f"Jet_eta_probe_{system}", f"B_{system}", "weight"),
                    ]
                )

        return self
    
    
    def __sample_cut(self, rdf : RNode) -> RNode:
        min_pt = 15
        tag_eta = 1.3
        asymmetry_alpha = 0.7
        delta_phi = 2.7
        rdf_dijet = (rdf.Filter("nJet >= 2")
                    # Choose tag and probe jets, use gRandom-Rndm()
                    .Define("tag_idx", "int(round(gRandom->Rndm()))")
                    .Define("probe_idx", "tag_idx == 0")

                    .Filter(f"Jet_pt[tag_idx] > {min_pt} && Jet_pt[probe_idx] > {min_pt}")
                    .Filter("Jet_jetId[tag_idx] >= 4 && Jet_jetId[probe_idx] >= 4")
                    .Define("deltaPhi_dijet", "abs(Jet_phi[tag_idx] - Jet_phi[probe_idx])")
                    .Redefine("deltaPhi_dijet", "deltaPhi_dijet > TMath::Pi() ? TMath::TwoPi() - deltaPhi_dijet : deltaPhi_dijet")
                    .Filter(f"abs(Jet_eta[tag_idx]) < {tag_eta} && deltaPhi_dijet > {delta_phi}")
                    # Vectors for the tag and probe jets
                    .Define("Dijet_tagVector", "ROOT::Math::PtEtaPhiMVector(Jet_pt[tag_idx], Jet_eta[tag_idx], Jet_phi[tag_idx], Jet_mass[tag_idx])")
                    .Define("Dijet_tagPolar", "ROOT::Math::Polar2DVectorF(Jet_pt[tag_idx], Jet_phi[tag_idx])")
                    .Define("Dijet_probePolar", "ROOT::Math::Polar2DVectorF(Jet_pt[probe_idx], Jet_phi[probe_idx])")
                    .Define("Dijet_bisectorPolar", "ROOT::Math::Polar2DVectorF(1.0, (Dijet_tagPolar - Dijet_probePolar).Phi())")
                    .Define("Dijet_probeVector", "ROOT::Math::PtEtaPhiMVector(Jet_pt[probe_idx], Jet_eta[probe_idx], Jet_phi[probe_idx], Jet_mass[probe_idx])")
                    .Define("Dijet_metPolar", "ROOT::Math::Polar2DVectorF(RawPuppiMET_pt, RawPuppiMET_phi)")
                    .Define("deltaEta_dijet", "abs(Jet_eta[tag_idx] - Jet_eta[probe_idx])")
                    .Define("deltaR_dijet", "sqrt(deltaPhi_dijet*deltaPhi_dijet + deltaEta_dijet*deltaEta_dijet)")
                    .Define("average_Pt_dijet", "(Jet_pt[tag_idx] + Jet_pt[probe_idx]) / 2.0")
                    # Filter the two jets
                    .Filter(f"(abs(Jet_pt[tag_idx] - Jet_pt[probe_idx]) / (2.0 * average_Pt_dijet)) < {asymmetry_alpha}")
                    .Define("Jet_pt_tag", "Jet_pt[tag_idx]")
                    .Define("Jet_pt_probe", "Jet_pt[probe_idx]")
                    .Define("Jet_eta_tag_dijet", "Jet_eta[tag_idx]") # Scuffed
                    .Define("Jet_eta_probe_dijet", "Jet_eta[probe_idx]")
                    .Define("Jet_phi_tag_dijet", "Jet_phi[tag_idx]")
                    .Define("Jet_phi_probe_dijet", "Jet_phi[probe_idx]")
                    # Asymmetry of the system
                    .Define("Asymmetry_dijet", "(Jet_pt[probe_idx] - Jet_pt[tag_idx]) / (Jet_pt[tag_idx] + Jet_pt[probe_idx])")
                    .Define("Dijet_ptAvp", "0.5 * (Dijet_tagPolar.Dot(Dijet_probePolar) - Dijet_probePolar.Dot(Dijet_probePolar))")
                    .Define("B_response", "1.0 + Dijet_metPolar.Dot(Dijet_bisectorPolar) / Dijet_ptAvp")
                    .Define("B_dijet", "RawPuppiMET_pt * cos(ROOT::VecOps::DeltaPhi(Jet_phi[tag_idx], RawPuppiMET_phi)) / (Jet_pt[tag_idx] + Jet_pt[probe_idx])")
                    # TODO: Add jet veto
                    )
        return rdf_dijet
    
        