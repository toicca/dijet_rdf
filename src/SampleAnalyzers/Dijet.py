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
 
    def do_sample_control(self) -> "DijetAnalyzer":
        print(f"Creating control histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            control_rdf = self.__sample_cut(self.Flag_cut(rdf))
            
            self.histograms[trigger].extend([
                control_rdf.Histo1D((f"Control_{self.system}_PtAvg", "Control_"+ str(self.system) + "_PtAvg;p_{T, ave} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_ptAvg", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_PtProbe", "Control_"+ str(self.system) + "_PtProbe;p_{T, probe} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_probePt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_PtTag", "Control_"+ str(self.system) + "_PtTag;p_{T, tag} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_tagPt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_EtaProbe", "Control_"+ str(self.system) + "_EtaProbe;#eta_{probe};N_{events}",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"]),
                                    "Dijet_probeEta", "weight"),
            ])
            
        return self
 

    def do_DB(self) -> "DijetAnalyzer":
        system = self.system
        print(f"Creating DB histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))    
            
            self.histograms[trigger].extend([
                # 3D Distributions from which information can be projected out
                db_rdf.Histo3D((f"DB_{system}_PtAvgVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, ave} (GeV);#eta_{probe};response;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_ptAvg", "Dijet_probeEta", "Dijet_dbResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtProbeVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, probe} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_probePt", "Dijet_probeEta", "Dijet_dbResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtTagVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, tag} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_tagPt", "Dijet_probeEta", "Dijet_dbResponse", "weight"),
                # 2D Asymmetry histogram for veto maps
                db_rdf.Profile2D((f"DB_{system}_EtaprobeVsPhiprobeVsAsymmetry", "DB_"+ str(system) + "_EtaprobeVsPhiprobeVsAsymmetry;#eta_{probe};#phi_{probe};asymmetry;N_{events}",
                                  self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                  "Dijet_probeEta", "Dijet_probePhi", "Dijet_A", "weight"),
            ])


        return self
    
    def do_MPF(self) -> "DijetAnalyzer":
        system = self.system
        print(f"Creating MPF histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            mpf_rdf = self.__sample_cut(self.Flag_cut(rdf))
            mpf_barrel_rdf = mpf_rdf.Filter("abs(Dijet_probeEta) < 1.3")        
            
            self.histograms[trigger].extend([
                # 3D Distributions
                mpf_rdf.Histo3D((f"MPF_{system}_PtAvgVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, ave} (GeV);#eta_{probe};response;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_ptAvg", "Dijet_probeEta", "Dijet_mpfResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtProbeVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, probe} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_probePt", "Dijet_probeEta", "Dijet_mpfResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtTagVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, tag} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_tagPt", "Dijet_probeEta", "Dijet_mpfResponse", "weight"),
                # Average distributions for pT in barrel
                # mpf_barrel_rdf.Profile1D((f"MPF_{system}_PtAveVsResponse", "MPF_"+ str(system) + "_PtAveVsResponse;p_{T, ave} (GeV);response;N_{events}",
                #                         self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                #                         "Dijet_ptAvg", "Dijet_mpfResponse", "weight"),
                # mpf_barrel_rdf.Profile1D((f"MPF_{system}_PtProbeVsResponse", "MPF_"+ str(system) + "_PtProbeVsResponse;p_{T, probe} (GeV);response;N_{events}",
                #                         self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                #                         "Dijet_probePt", "Dijet_mpfResponse", "weight"),
                # mpf_barrel_rdf.Profile1D((f"MPF_{system}_PtTagVsResponse", "MPF_"+ str(system) + "_PtTagVsResponse;p_{T, tag} (GeV);response;N_{events}",
                #                         self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                #                         "Dijet_tagPt", "Dijet_mpfResponse", "weight"),
            ])

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
                    .Define("Dijet_ptAvg", "(Jet_pt[tag_idx] + Jet_pt[probe_idx]) / 2.0")
                    # Filter the two jets
                    .Filter(f"(abs(Jet_pt[tag_idx] - Jet_pt[probe_idx]) / (2.0 * Dijet_ptAvg)) < {asymmetry_alpha}")
                    .Define("Dijet_tagPt", "Jet_pt[tag_idx]")
                    .Define("Dijet_probePt", "Jet_pt[probe_idx]")
                    .Define("Dijet_tagEta", "Jet_eta[tag_idx]") # Scuffed
                    .Define("Dijet_probeEta", "Jet_eta[probe_idx]")
                    .Define("Dijet_tagPhi", "Jet_phi[tag_idx]")
                    .Define("Dijet_probePhi", "Jet_phi[probe_idx]")
                    # Responses
                    .Define("Dijet_A", "(Jet_pt[probe_idx] - Jet_pt[tag_idx]) / (2.0 * Dijet_ptAvg)")
                    .Define("Dijet_B", "Dijet_metPolar.Dot(Dijet_tagPolar) / (2.0 * Dijet_ptAvg * Dijet_tagPolar.R())")
                    .Define("Dijet_ptAvp", "0.5 * (Dijet_tagPolar.Dot(Dijet_probePolar) - Dijet_probePolar.Dot(Dijet_probePolar))")
                    .Define("Dijet_mpfResponse", "(1.0 + Dijet_B) / (1.0 - Dijet_B)")
                    .Define("Dijet_dbResponse", "(1.0 + Dijet_A) / (1.0 - Dijet_A)")
                    )
        return rdf_dijet
    
        