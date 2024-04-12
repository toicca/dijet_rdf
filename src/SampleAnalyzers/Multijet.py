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
    
    def do_sample_control(self) -> "MultijetAnalyzer":
        print(f"Creating control histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            control_rdf = self.__sample_cut(self.Flag_cut(rdf))
            
            self.histograms[trigger].extend([
                control_rdf.Histo1D((f"Control_{self.system}_PtAvg", "Control_"+ str(self.system) + "_PtAvg;p_{T, ave} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_ptAvg", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_PtRecoil", "Control_"+ str(self.system) + "_PtRecoil;p_{T, recoil} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_recoilPt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_PtLead", "Control_"+ str(self.system) + "_PtLead;p_{T, lead} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_leadPt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_EtaRecoil", "Control_"+ str(self.system) + "_EtaRecoil;#eta_{Recoil};N_{events}",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"]),
                                    "Multijet_recoilEta", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_nRecoil", "Control_"+ str(self.system) + "_nRecoil;N_{recoil};N_{events}",
                                    20, 0, 20),
                                    "Multijet_nRecoil", "weight"),
            ])
            
        return self

    def do_DB(self) -> "MultijetAnalyzer":
        system = "multijet"
        print(f"Creating DB histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))
            db_barrel_rdf = db_rdf.Filter("abs(Multijet_recoilEta) < 1.3")        
            
            self.histograms[trigger].extend([
                # 3D Distributions
                db_rdf.Histo3D((f"DB_{system}_PtAvgVsEtaVsResponse", "DB_"+ str(system) + "_PtVsEtaVsResponse;p_{T, ave} (GeV);#eta_{recoil};response;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_ptAvg", "Multijet_recoilEta", "Multijet_dbResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtRecoilVsEtaVsResponse", "DB_"+ str(system) + "_PtVsEtaVsResponse;p_{T, lead} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_dbResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtLeadVsEtaVsResponse", "DB_"+ str(system) + "_PtVsEtaVsResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_dbResponse", "weight"),
                # Average distributions for pT in barrel
                db_barrel_rdf.Profile1D((f"DB_{system}_PtAveVsResponse", "DB_"+ str(system) + "_PtAveVsResponse;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvg", "Multijet_dbResponse", "weight"),
                db_barrel_rdf.Profile1D((f"DB_{system}_PtLeadVsResponse", "DB_"+ str(system) + "_PtRecoilVsResponse;p_{T, lead} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_dbResponse", "weight"),
                db_barrel_rdf.Profile1D((f"DB_{system}_PtRecoilVsResponse", "DB_"+ str(system) + "_PtLeadVsResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_leadPt", "Multijet_dbResponse", "weight"),
                db_barrel_rdf.Profile1D((f"DB_{system}_PtAveVsR", "DB_"+ str(system) + "_PtAveVsResponse;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvg", "Multijet_dbR", "weight"),
                db_barrel_rdf.Profile1D((f"DB_{system}_PtRecoilVsR", "DB_"+ str(system) + "_PtRecoilVsResponse;p_{T, lead} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_dbR", "weight"),
                db_barrel_rdf.Profile1D((f"DB_{system}_PtLeadVsR", "DB_"+ str(system) + "_PtLeadVsResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_leadPt", "Multijet_dbR", "weight"),
            ])

        return self
    
    def do_MPF(self) -> "MultijetAnalyzer":
        system = "multijet"
        print(f"Creating MPF histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            mpf_rdf = self.__sample_cut(self.Flag_cut(rdf))
            mpf_barrel_rdf = mpf_rdf#.Filter("abs(Multijet_recoilEta) < 1.3")        
            
            self.histograms[trigger].extend([
                # 3D Distributions
                mpf_rdf.Histo3D((f"MPF_{system}_PtAvgVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, ave} (GeV);#eta_{recoil};response;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_ptAvg", "Multijet_recoilEta", "Multijet_mpfResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtRecoilVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, lead} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_mpfResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtLeadVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_mpfResponse", "weight"),
                # Average distributions for pT in barrel
                mpf_barrel_rdf.Profile1D((f"MPF_{system}_PtAveVsResponse", "MPF_"+ str(system) + "_PtAveVsResponse;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvg", "Multijet_mpfAvgResponse", "weight"),
                mpf_barrel_rdf.Profile1D((f"MPF_{system}_PtAvpVsResponse", "MPF_"+ str(system) + "_PtAvpVsResponse;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvp", "Multijet_mpfAvpResponse", "weight"),
                mpf_barrel_rdf.Profile1D((f"MPF_{system}_PtLeadVsResponse", "MPF_"+ str(system) + "_PtLeadVsResponse;p_{T, lead} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_leadPt", "Multijet_mpfLeadResponse", "weight"),
                mpf_barrel_rdf.Profile1D((f"MPF_{system}_PtRecoilVsResponse", "MPF_"+ str(system) + "_PtRecoilVsResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_mpfRecoilResponse", "weight"),
            ])

        return self
    
    def __sample_cut(self, rdf : RNode) -> RNode:
        min_pt = 30
        lead_eta = 1.3
        recoil_eta = 2.5
        asymmetry_alpha = 0.7
        delta_phi = 1.0
        
        if not hasattr(ROOT, "sum_as_four_vectors"):
            ROOT.gInterpreter.Declare("""
                                    ROOT::Math::PtEtaPhiMVector sum_as_four_vectors(ROOT::RVec<float> pt,ROOT:: RVec<float> eta, ROOT::RVec<float> phi, ROOT::RVec<float> mass) {
                                        ROOT::Math::PtEtaPhiMVector sum(0, 0, 0, 0);
                                        for (int i = 0; i < pt.size(); i++) {
                                            sum += ROOT::Math::PtEtaPhiMVector(pt[i], eta[i], phi[i], mass[i]);
                                        }
                                        return sum;
                                    }
                                    """)
                                  
        
        rdf_multijet = (rdf.Filter("nJet >= 3")
                    # .Filter(f"Jet_pt[Jet_order[0]] > {min_pt} && abs(Jet_eta[Jet_order[0]]) < {lead_eta} && Jet_jetId[Jet_order[0]] >= 4")
                    .Define("Jet_eta_new", f"Jet_eta[Jet_pt > {min_pt} && abs(Jet_eta) < {recoil_eta} && Jet_jetId >= 4]")
                    .Redefine("Jet_phi", f"Jet_phi[Jet_pt > {min_pt} && abs(Jet_eta) < {recoil_eta} && Jet_jetId >= 4]")
                    .Redefine("Jet_mass", f"Jet_mass[Jet_pt > {min_pt} && abs(Jet_eta) < {recoil_eta} && Jet_jetId >= 4]")
                    .Define("Jet_jetIdNew", f"Jet_jetId[Jet_pt > {min_pt} && abs(Jet_eta) < {recoil_eta} && Jet_jetId >= 4]")
                    .Redefine("Jet_order", f"Jet_order[Jet_pt > {min_pt} && abs(Jet_eta) < {recoil_eta} && Jet_jetId >= 4]")
                    .Redefine("Jet_pt", f"Jet_pt[Jet_pt > {min_pt} && abs(Jet_eta) < {recoil_eta} && Jet_jetId >= 4]")
                    .Redefine("Jet_eta", "Jet_eta_new")
                    .Redefine("Jet_jetId", "Jet_jetIdNew")
                    .Filter("Jet_pt.size() >= 3")
                    .Define("leading_idx", "Jet_order[0]")
                    .Define("second_idx", "Jet_order[1]")
                    .Define("third_idx", "Jet_order[2]")
                    .Define("Multijet_leadPt", "Jet_pt[Jet_order[0]]")
                    .Define("Jet_pt_recoil", f"Jet_pt[Jet_order > 0 && (ROOT::VecOps::DeltaPhi(ROOT::RVec<float>(Jet_phi), float(Jet_phi[Jet_order[0]])) > {delta_phi})]")
                    .Define("Jet_eta_recoil", f"Jet_eta[Jet_order > 0 && (ROOT::VecOps::DeltaPhi(ROOT::RVec<float>(Jet_phi), float(Jet_phi[Jet_order[0]])) > {delta_phi})]")
                    .Define("Jet_mass_recoil", f"Jet_mass[Jet_order > 0 && (ROOT::VecOps::DeltaPhi(ROOT::RVec<float>(Jet_phi), float(Jet_phi[Jet_order[0]])) > {delta_phi})]")
                    .Define("Jet_phiNew_recoil", f"Jet_phi[Jet_order > 0 && (ROOT::VecOps::DeltaPhi(ROOT::RVec<float>(Jet_phi), float(Jet_phi[Jet_order[0]])) > {delta_phi})]")
                    .Define("Jet_order_recoil", f"Jet_order[Jet_order > 0 && (ROOT::VecOps::DeltaPhi(ROOT::RVec<float>(Jet_phi), float(Jet_phi[Jet_order[0]])) > {delta_phi})]")
                    .Define("Jet_phi_recoil", "Jet_phiNew_recoil")
                    .Define("Multijet_recoilVector", "sum_as_four_vectors(Jet_pt_recoil, Jet_eta_recoil, Jet_phi_recoil, Jet_mass_recoil)")
                    .Define("Multijet_recoilPt", "Multijet_recoilVector.Pt()")
                    .Define("Multijet_recoilEta", "Multijet_recoilVector.Eta()")
                    .Define("Multijet_recoilPhi", "Multijet_recoilVector.Phi()")
                    .Define("Multijet_recoilMass", "Multijet_recoilVector.M()")
                    .Redefine("second_idx", "Jet_order_recoil[0]")
                    .Redefine("third_idx", "Jet_order_recoil[1]")
                    .Filter(f"Jet_pt[second_idx] < 0.6 * Multijet_recoilPt")
                    .Filter(f"Jet_pt[third_idx] < 0.6 * Multijet_recoilPt")
                    # Vectors
                    .Define("Multijet_leadVector", "ROOT::Math::PtEtaPhiMVector(Jet_pt[leading_idx], Jet_eta[leading_idx], Jet_phi[leading_idx], Jet_mass[leading_idx])")
                    .Define("Multijet_leadPolar", "ROOT::Math::Polar2DVectorF(Jet_pt[leading_idx], Jet_phi[leading_idx])")
                    .Define("Multijet_recoilPolar", "ROOT::Math::Polar2DVectorF(Multijet_recoilPt, Multijet_recoilPhi)")
                    .Define("Multijet_avgPolar", "(Multijet_recoilPolar- Multijet_leadPolar) / (Multijet_leadPolar.R() + Multijet_recoilPolar.R())")
                    .Define("Multijet_bisectorPolar", "ROOT::Math::Polar2DVectorF(1.0, (Multijet_recoilPolar - Multijet_leadPolar).Phi())")
                    .Define("Multijet_nRecoil", "Jet_pt.size() - 1")
                    .Define("Multijet_ptAvg", "(Multijet_leadPt + Multijet_recoilPt) / float(Multijet_nRecoil + 1)")
                    .Define("Multijet_ptAvp", "0.5 * (Multijet_recoilPolar.Dot(Multijet_bisectorPolar) - Multijet_leadPolar.Dot(Multijet_bisectorPolar))")
                    # Responses
                    .Define("Multijet_A", "(Multijet_leadPt - Multijet_recoilPt) / (Multijet_leadPt + Multijet_recoilPt)")
                    .Define("Multijet_B", "RawPuppiMET_polar.Dot(Multijet_leadPolar) / (2.0 * Multijet_ptAvg * Multijet_leadPolar.R())")
                    .Define("Multijet_mpfResponse", "(1.0 + Multijet_B) / (1.0 - Multijet_B)")
                    .Define("Multijet_dbResponse", "(1.0 + Multijet_A) / (1.0 - Multijet_A)")
                    .Define("Multijet_dbR", "Multijet_leadPolar.R() / Multijet_recoilPolar.R()")
                    .Define("Multijet_mpfRecoilResponse", "1.0 + RawPuppiMET_polar.Dot(Multijet_recoilPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")
                    .Define("Multijet_mpfLeadResponse", "1.0 + RawPuppiMET_polar.Dot(Multijet_leadPolar) / (Multijet_leadPolar.R() * Multijet_leadPolar.R())")
                    .Define("Multijet_mpfAvgResponse", "1.0 + RawPuppiMET_polar.Dot(Multijet_avgPolar) / (Multijet_ptAvg * Multijet_avgPolar.R())")
                    .Define("Multijet_mpfAvpResponse", "1.0 + RawPuppiMET_polar.Dot(Multijet_bisectorPolar) / (Multijet_ptAvg)")

                    .Define("Multijet_deltaPhi", "ROOT::VecOps::DeltaPhi(Multijet_leadVector.Phi(), Multijet_recoilVector.Phi())")
                    # .Filter(f"abs(Multijet_deltaPhi - TMath::Pi()) < 0.3")
                    
                    )

        return rdf_multijet
        