import ROOT
from typing import List
import numpy as np
from RDFAnalyzer import RDFAnalyzer, JEC_corrections
    
RDataFrame = ROOT.RDataFrame
RunGraphs = ROOT.RDF.RunGraphs
RNode = ROOT.RDF.RNode
    
class MultijetAnalyzer(RDFAnalyzer):
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
        super().__init__(filelist, trigger_list, json_file, nFiles, JEC, nThreads, progress_bar, isMC=isMC, local=local, system="multijet")
 
    def Flag_cut(self, rdf: RNode) -> RNode:
        return super().Flag_cut(rdf)

    def do_DB(self) -> "MultijetAnalyzer":
        system = "multijet"
        print(f"Creating DB histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))
            pT_binLabels = ["average_Pt_multijet", "Jet_pt_lead", "pt_recoil"]

            db_rdf = (db_rdf.Define(f"response_DB_{system}", f"(1.0 + Asymmetry_{system}) / (1.0 - Asymmetry_{system})")
                    )
            
            
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"DB_{system}_Response", "DB_"+ str(system) + "_response;response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), f"response_DB_{system}", "weight"),
                db_rdf.Histo2D((f"DB_{system}_EtaprobeVsResponse", "DB_"+ str(system) + "_EtaVsResponse;|#eta_{probe}|;response",
                                self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                               f"pt_recoil", f"response_DB_{system}", "weight"),
                db_rdf.Histo1D((f"DB_{system}_Asymmetry", "DB_"+ str(system) + "_Asymmetry;Asymmetry;N_{events}", self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]), 
                               f"Asymmetry_{system}", "weight"),
                db_rdf.Profile2D((f"DB_{system}_EtaprobeVsPhiprobeVsAsymmetry", "DB_"+ str(system) + "_EtaVsPhiVsAsymmetry;|#eta_{probe}|;#phi_{probe};Asymmetry", 
                                self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]), 
                                f"pt_recoil", f"phi_recoil", f"Asymmetry_{system}", "weight"),
            ])
            for label in pT_binLabels:
                self.histograms[trigger].append(
                    db_rdf.Histo2D((f"DB_{system}_PtVsResponse_{label}", "DB_"+str(system)+"_PtVsResponse_"+label+";p_{T, "+label+"} (GeV);response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["response"]["n"], self.bins["response"]["bins"]),
                                    label, f"response_DB_{system}", "weight"),
                )
        return self
    
    def do_MPF(self) -> "MultijetAnalyzer":
        system = "multijet"
        print(f"Creating MPF histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))
            pT_binLabels = ["average_Pt_multijet", "Jet_pt_lead", "pt_recoil"]

            db_rdf = (db_rdf.Define(f"response_MPF_{system}", f"(1.0 + B_{system}) / (1.0 - B_{system})")
                    )
            
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"MPF_{system}_Response", "MPF_"+ str(system) + "_response;response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), f"response_MPF_{system}", "weight"),
                db_rdf.Histo2D((f"DB_{system}_EtaVsResponse", "DB_"+ str(system) + "_EtaVsResponse;|#eta|;response", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), f"Jet_eta_tag_{system}", f"response_MPF_{system}", "weight")
            ])
            for label in pT_binLabels:
                self.histograms[trigger].append(
                    db_rdf.Histo2D((f"MPF_{system}_PtVsResponse_{label}", "MPF_"+str(system)+"_PtVsResponse_"+label+";p_{T, "+label+"} (GeV);response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["response"]["n"], self.bins["response"]["bins"]),
                                    label, f"response_MPF_{system}", "weight"),
                )

        return self
    
    def __sample_cut(self, rdf : RNode) -> RNode:
        min_pt = 30
        lead_eta = 1.3
        recoil_eta = 2.5
        asymmetry_alpha = 0.7
        delta_phi = 0.3
        rdf_dijet = (rdf.Filter("nJet >= 3")
                    .Redefine("Jet_eta", f"Jet_eta[Jet_pt > {min_pt}]")
                    .Redefine("Jet_phi", f"Jet_phi[Jet_pt > {min_pt}]")
                    .Redefine("Jet_mass", f"Jet_mass[Jet_pt > {min_pt}]")
                    .Redefine("Jet_jetId", f"Jet_jetId[Jet_pt > {min_pt}]")
                    .Redefine("Jet_pt", f"Jet_pt[Jet_pt > {min_pt}]")
                    .Filter("Jet_pt.size() >= 3")
                    .Define("leading_idx", "Jet_order[0]")
                    .Define("second_idx", "Jet_order[1]")
                    .Define("third_idx", "Jet_order[2]")
                    # Filter the lead and recoil system
                    .Define("Jet_pt_lead", "Jet_pt[leading_idx]")
                    .Filter(f"Jet_pt_lead > {min_pt} && abs(Jet_eta[leading_idx]) < {lead_eta} && Jet_jetId[leading_idx] >= 4")
                    .Filter(f"Jet_pt[second_idx] > {min_pt} && abs(Jet_eta[second_idx]) < {recoil_eta} && Jet_jetId[second_idx] >= 4")
                    .Filter(f"Jet_pt[third_idx] > {min_pt} && abs(Jet_eta[third_idx]) < {recoil_eta} && Jet_jetId[third_idx] >= 4")
                    
                    .Define("pt_recoil", "ROOT::VecOps::Sum(Jet_pt) - Jet_pt_lead")
                    .Define("eta_recoil", "ROOT::VecOps::Sum(Jet_eta) - Jet_eta[leading_idx]")
                    .Define("phi_recoil", "ROOT::VecOps::Sum(Jet_phi) - Jet_phi[leading_idx]")
                    .Define("mass_recoil", "ROOT::VecOps::Sum(Jet_mass) - Jet_mass[leading_idx]")
                    .Define("average_Pt_multijet", "(Jet_pt_lead + pt_recoil) / (float(Jet_pt.size()))")
                    
                    .Define("Jet_eta_tag_multijet", "eta_recoil") # Also, scuffed
                    .Define("deltaPhi_multijet", "abs(Jet_phi[leading_idx] - phi_recoil)")
                    .Redefine("deltaPhi_multijet", "deltaPhi_multijet > TMath::Pi() ? TMath::TwoPi() - deltaPhi_multijet : deltaPhi_multijet")

                    .Filter(f"Jet_pt[second_idx] < 0.6 * pt_recoil")
                    .Filter(f"Jet_pt[third_idx] < 0.6 * pt_recoil")
                    .Filter(f"abs(deltaPhi_multijet - TMath::Pi()) < {delta_phi}")
                    # Asymmetry of the system
                    .Define("Asymmetry_multijet", "(Jet_pt_lead - pt_recoil) / (Jet_pt_lead + pt_recoil)")
                    .Define("B_multijet", "RawPuppiMET_pt * cos(ROOT::VecOps::DeltaPhi(Jet_phi[leading_idx], RawPuppiMET_phi)) / (Jet_pt_lead + pt_recoil)")
                    # TODO: Add jet veto
                    )
        # self.histograms["all"].extend([rdf_dijet.Histo1D(("Asymmetry", "Asymmetry;Asymmetry;N_{events}", 100, -1.0, 1.0), f"Asymmetry_multijet", "weight"),
        #                                rdf_dijet.Histo1D(("pt_recoil", "pt_recoil;p_{T, recoil} (GeV);N_{events}", 100, 0.0, 1000.0), "pt_recoil", "weight"),
        #                                rdf_dijet.Histo1D(("Jet_pt_lead", "Jet_pt_lead;p_{T, lead} (GeV);N_{events}", 100, 0.0, 1000.0), "Jet_pt_lead", "weight"),
        #                                rdf_dijet.Histo1D(("deltaPhi_multijet", "deltaPhi_multijet;#Delta#phi_{lead, recoil};N_{events}", 100, 0.0, 3.141592653), "deltaPhi_multijet", "weight"),
        #                                rdf_dijet.Histo1D(("average_Pt_multijet", "average_Pt_multijet;average p_{T} (GeV);N_{events}", 100, 0.0, 1000.0), "average_Pt_multijet", "weight"),])

        return rdf_dijet
        