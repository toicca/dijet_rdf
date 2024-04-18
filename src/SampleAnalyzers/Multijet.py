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
                local : bool = False,
                run_raw : bool = False,
                selection_only : bool = False
                ):
        super().__init__(filelist, trigger_list, json_file, nFiles, JEC, nThreads, progress_bar, isMC=isMC, local=local, system="multijet", run_raw=run_raw, selection_only=selection_only)
 
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
                control_rdf.Histo1D((f"Control_{self.system}_leadIdx", "Control_"+ str(self.system) + "_leadIdx;leadIdx;N_{events}",
                                    20, 0, 20),
                                    "lead_idx", "weight"),
                
                # MET direction as a function of pT
                control_rdf.Profile1D((f"Control_{self.system}_METdotRecoilVsPtRecoil", "Control_"+ str(self.system) + "_METdotRecoilVsPtRecoil;p_{T, recoil} (GeV);MET dot recoil;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_recoilPt", "Multijet_METdotRecoil", "weight"),
                control_rdf.Profile1D((f"Control_{self.system}_METdotLeadVsPtLead", "Control_"+ str(self.system) + "_METdotLeadVsPtLead;p_{T, lead} (GeV);MET dot lead;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_leadPt", "Multijet_METdotLead", "weight"),
                
                # Deltaphi to lead jet
                control_rdf.Histo1D((f"Control_{self.system}_DeltaPhi", "Control_"+ str(self.system) + "_DeltaPhi;#Delta#phi;N_{events}",
                                    20, 0, 3.2),
                                    "Multijet_deltaPhi", "weight"),
                
                # Multijet asymmetry map
                control_rdf.Profile2D((f"Control_{self.system}_EtaRVsPhiRVsAsymmetry", "Control_"+ str(self.system) + "_EtaRVsPhiRVsAsymmetry;#eta_{recoil};#phi_{recoil};asymmetry;N_{events}",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                    "Multijet_recoilEta", "Multijet_recoilPhi", "Multijet_A", "weight"),
                
                control_rdf.Profile1D((f"Control_{self.system}_PtRecoilVsNRecoil", "MPF_"+ str(self.system) + "_PtRecoilVsNRecoil;p_{T, recoil} (GeV);N_{recoil};N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_nRecoil", "weight"),
            ])
            
        return self

    def do_DB(self) -> "MultijetAnalyzer":
        system = "multijet"
        print(f"Creating DB histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))    
            
            self.histograms[trigger].extend([
                # 3D Distributions
                db_rdf.Histo3D((f"DB_{system}_PtAvgVsEtaVsResponse", "DB_"+ str(system) + "_PtVsEtaVsResponse;p_{T, ave} (GeV);#eta_{recoil};response;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_ptAvg", "Multijet_recoilEta", "Multijet_dbResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtRecoilVsEtaVsResponse", "DB_"+ str(system) + "_PtVsEtaVsResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_dbResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtLeadVsEtaVsResponse", "DB_"+ str(system) + "_PtVsEtaVsResponse;p_{T, lead} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_dbResponse", "weight"),
                # 3D Distributions for A
                db_rdf.Histo3D((f"DB_{system}_PtAvgVsEtaVsA", "DB_"+ str(system) + "_PtAvgVsEtaVsA;p_{T, ave} (GeV);#eta_{recoil};asymmetry;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Multijet_ptAvg", "Multijet_recoilEta", "Multijet_A", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtRecoilVsEtaVsA", "DB_"+ str(system) + "_PtRecoilVsEtaVsA;p_{T, recoil} (GeV);#eta_{recoil};asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_A", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtLeadVsEtaVsA", "DB_"+ str(system) + "_PtLeadVsEtaVsA;p_{T, lead} (GeV);#eta_{recoil};asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_A", "weight"),
                
                
                # Average distributions for pT in barrel
                db_rdf.Profile1D((f"DB_{system}_PtAveVsResponse", "DB_"+ str(system) + "_PtAveVsResponse;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvg", "Multijet_dbResponse", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtRecoilVsResponse", "DB_"+ str(system) + "_PtRecoilVsResponse;p_{T, lead} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_dbResponse", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtLeadVsResponse", "DB_"+ str(system) + "_PtLeadVsResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_leadPt", "Multijet_dbResponse", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtAveVsR", "DB_"+ str(system) + "_PtAveVsR;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvg", "Multijet_dbR_avg", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtRecoilVsR_b2b", "DB_"+ str(system) + "_PtRecoilVsR_b2b;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_dbR_recoilB2B", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtRecoilVsR", "DB_"+ str(system) + "_PtRecoilVsR;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_dbR_recoil", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtLeadVsR", "DB_"+ str(system) + "_PtLeadVsR;p_{T, lead} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_leadPt", "Multijet_dbR_lead", "weight"),
            ])
            db_gt100_rdf = self.__sample_cut(self.Flag_cut(rdf)).Filter("Multijet_recoilPt > 100")
            self.histograms[trigger].extend([
                db_gt100_rdf.Profile1D((f"DB_{system}_RunVsResponse", "DB_"+ str(system) + "_RunVsResponse;Run;response;N_{events}",
                                         self.bins["runs"]["n"], self.bins["runs"]["bins"]),
                                            "run", "Multijet_dbResponse", "weight"),
            ])
        return self
    
    def do_MPF(self) -> "MultijetAnalyzer":
        system = "multijet"
        print(f"Creating MPF histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            mpf_rdf = self.__sample_cut(self.Flag_cut(rdf))   
            
            self.histograms[trigger].extend([
                # 3D Distributions for response
                mpf_rdf.Histo3D((f"MPF_{system}_PtAvgVsEtaVsResponse", "MPF_"+ str(system) + "_PtAvgVsEtaVsResponse;p_{T, ave} (GeV);#eta_{recoil};response;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_ptAvg", "Multijet_recoilEta", "Multijet_mpfAvgResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtRecoilVsEtaVsResponse", "MPF_"+ str(system) + "_PtRecoilVsEtaVsResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_mpfRecoilResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtLeadVsEtaVsResponse", "MPF_"+ str(system) + "_PtLeadVsEtaVsResponse;p_{T, lead} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_mpfLeadResponse", "weight"),
                # 3D Distributions for B
                mpf_rdf.Histo3D((f"MPF_{system}_PtAvgVsEtaVsB", "MPF_"+ str(system) + "_PtAvgVsEtaVsB;p_{T, ave} (GeV);#eta_{recoil};response;N_{events}",
                                 self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Multijet_ptAvg", "Multijet_recoilEta", "Multijet_B", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtRecoilVsEtaVsB", "MPF_"+ str(system) + "_PtRecoilVsEtaVsB;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_Brecoil", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtLeadVsEtaVsB", "MPF_"+ str(system) + "_PtLeadVsEtaVsB;p_{T, lead} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_B", "weight"),

                # Average distributions for pT in barrel
                mpf_rdf.Profile1D((f"MPF_{system}_PtAveVsResponse", "MPF_"+ str(system) + "_PtAveVsResponse;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvg", "Multijet_mpfAvgResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtAvpVsResponse", "MPF_"+ str(system) + "_PtAvpVsResponse;p_{T, ave} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_ptAvp", "Multijet_mpfAvpResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtLeadVsResponse", "MPF_"+ str(system) + "_PtLeadVsResponse;p_{T, lead} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_leadPt", "Multijet_mpfLeadResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtRecoilVsResponse", "MPF_"+ str(system) + "_PtRecoilVsResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_mpfRecoilResponse", "weight"),
            ])
            mpf_gt100_rdf = self.__sample_cut(self.Flag_cut(rdf)).Filter("Multijet_recoilPt > 100")
            self.histograms[trigger].extend([
                mpf_gt100_rdf.Profile1D((f"MPF_{system}_RunVsResponse", "MPF_"+ str(system) + "_RunVsResponse;Run;response;N_{events}",
                                         self.bins["runs"]["n"], self.bins["runs"]["bins"]),
                                            "run", "Multijet_mpfRecoilResponse", "weight"),
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
            
            ROOT::RVec<float> DeltaR_to_lead(float lead_eta, ROOT::RVec<float> recoil_eta, float lead_phi, ROOT::RVec<float> recoil_phi){
                ROOT::RVec<float> result(recoil_eta.size());
                
                for (int i; i < result.size(); i++){
                    result[i] = ROOT::VecOps::DeltaR(lead_eta, recoil_eta[i], lead_phi, recoil_phi[i]);
                }
                
                return result;
            }
            
            ROOT::RVec<float> DeltaPhi_to_lead(float lead_phi, ROOT::RVec<float> recoil_phi){
                ROOT::RVec<float> result(recoil_phi.size());
                
                for (int i; i < result.size(); i++){
                    result[i] = ROOT::VecOps::DeltaPhi(lead_phi, recoil_phi[i]);
                }
                
                return result;
            }
                                    """)
                                  
        # TODO: p4l1rc vector for each recoil to do p4l1rc - p4 MET
        rdf_multijet = (rdf.Filter("nJet >= 3")
                    .Filter(f"Jet_pt[Jet_order[0]] > {min_pt} && abs(Jet_eta[Jet_order[0]]) < {lead_eta} && Jet_jetId[Jet_order[0]] >= 4")
                    .Filter(f"Jet_pt[Jet_order[1]] > {min_pt} && abs(Jet_eta[Jet_order[1]]) < {recoil_eta} && Jet_jetId[Jet_order[1]] >= 4")
                    .Filter(f"Jet_pt[Jet_order[2]] > {min_pt} && abs(Jet_eta[Jet_order[2]]) < {recoil_eta} && Jet_jetId[Jet_order[2]] >= 4")
                    .Define("lead_idx", "Jet_order[0]")
                    .Define("Multijet_leadPt", "Jet_pt[Jet_order[0]]")
                    .Define("Multijet_deltaPhi_to_lead", "DeltaPhi_to_lead(Jet_phi[lead_idx], Jet_phi)")
                    .Filter(f"ROOT::VecOps::All(abs(Jet_eta[(Jet_order != lead_idx) && (Jet_pt >= {min_pt})]) <= {recoil_eta})") # events without forward jets
                    .Filter(f"ROOT::VecOps::All(Multijet_deltaPhi_to_lead[(Jet_order != lead_idx) && (Jet_pt >= {min_pt})] >= {delta_phi})") # events with jets back-to-back 
                    .Define("Jet_pt_recoil", f"Jet_pt[(Jet_order != lead_idx) && (Jet_pt >= {min_pt}) && (abs(Jet_eta) <= {recoil_eta}) && (Multijet_deltaPhi_to_lead > 1.0)]")
                    .Define("Jet_eta_recoil", f"Jet_eta[(Jet_order != lead_idx) && (Jet_pt >= {min_pt}) && (abs(Jet_eta) <= {recoil_eta}) && (Multijet_deltaPhi_to_lead > 1.0)]")
                    .Define("Jet_mass_recoil", f"Jet_mass[(Jet_order != lead_idx) && (Jet_pt >= {min_pt}) && (abs(Jet_eta) <= {recoil_eta}) && (Multijet_deltaPhi_to_lead > 1.0)]")
                    .Define("Jet_phi_recoil", f"Jet_phi[(Jet_order != lead_idx) && (Jet_pt >= {min_pt}) && (abs(Jet_eta) <= {recoil_eta}) && (Multijet_deltaPhi_to_lead > 1.0)]")
                    .Define("Jet_order_recoil", f"Jet_order[(Jet_order != lead_idx) && (Jet_pt >= {min_pt}) && (abs(Jet_eta) <= {recoil_eta}) && (Multijet_deltaPhi_to_lead > 1.0)]")
                    .Filter("Jet_pt_recoil.size() >= 2")
                    .Define("Multijet_recoilVector", "sum_as_four_vectors(Jet_pt_recoil, Jet_eta_recoil, Jet_phi_recoil, Jet_mass_recoil)")
                    .Define("Multijet_recoilPt", "Multijet_recoilVector.Pt()")
                    .Filter(f"Jet_pt[Jet_order_recoil[0]] < 0.6 * Multijet_recoilPt")
                    .Filter(f"Jet_pt[Jet_order_recoil[1]] < 0.6 * Multijet_recoilPt")
                    .Define("Multijet_recoilEta", "Multijet_recoilVector.Eta()")
                    .Define("Multijet_recoilPhi", "Multijet_recoilVector.Phi()")
                    .Define("Multijet_recoilMass", "Multijet_recoilVector.M()")
                    # Vectors
                    .Define("Multijet_leadVector", "ROOT::Math::PtEtaPhiMVector(Jet_pt[lead_idx], Jet_eta[lead_idx], Jet_phi[lead_idx], Jet_mass[lead_idx])")
                    .Define("Multijet_leadPolar", "ROOT::Math::Polar2DVectorF(Jet_pt[lead_idx], Jet_phi[lead_idx])")
                    .Define("Multijet_recoilPolar", "ROOT::Math::Polar2DVectorF(Multijet_recoilPt, Multijet_recoilPhi)")
                    .Define("Multijet_avgPolar", "(Multijet_recoilPolar - Multijet_leadPolar) / 2.0")
                    .Define("Multijet_bisectorPolar", "ROOT::Math::Polar2DVectorF(1.0, (Multijet_recoilPolar + Multijet_leadPolar).Phi())")
                    .Define("Multijet_nRecoil", "Jet_pt_recoil.size()")
                    .Define("Multijet_ptAvg", "(Multijet_leadVector.Pt() + Multijet_recoilVector.Pt()) / 2.0")
                    .Define("Multijet_ptAvp", "0.5 * (Multijet_recoilPolar.Dot(Multijet_bisectorPolar) - Multijet_leadPolar.Dot(Multijet_bisectorPolar))")
                    .Define("Multijet_deltaPhi", "ROOT::VecOps::DeltaPhi(Multijet_leadVector.Phi(), Multijet_recoilVector.Phi())")
                    .Filter(f"abs(Multijet_deltaPhi - TMath::Pi()) < 0.3")
                    .Define("Multijet_METdotRecoil", "RawPuppiMET_polar.Dot(Multijet_recoilPolar) / (Multijet_recoilPolar.R() * RawPuppiMET_polar.R())")
                    .Define("Multijet_METdotLead", "RawPuppiMET_polar.Dot(Multijet_leadPolar) / (Multijet_leadPolar.R() * RawPuppiMET_polar.R())")
                    # Responses
                    # DB
                    .Define("Multijet_A", "(Multijet_recoilPt - Multijet_leadPt) / (Multijet_leadPt + Multijet_recoilPt)")
                    .Define("Multijet_dbResponse", "(1.0 + Multijet_A) / (1.0 - Multijet_A)")
                    .Define("Multijet_dbR_recoil", "-Multijet_recoilVector.Vect().Dot(Multijet_leadVector.Vect()) / (Multijet_recoilVector.Pt() * Multijet_recoilVector.Pt())") 
                    .Define("Multijet_dbR_recoilB2B", "Multijet_recoilVector.Pt() / Multijet_leadVector.Pt()")
                    .Define("Multijet_dbR_lead", "-Multijet_recoilVector.Vect().Dot(Multijet_leadVector.Vect()) / (Multijet_recoilVector.Pt() * Multijet_leadVector.Pt())")
                    .Define("Multijet_dbR_avg", "-Multijet_recoilVector.Vect().Dot(Multijet_leadVector.Vect()) / (Multijet_recoilVector.Pt() * Multijet_ptAvg)")
                    # MPF
                    .Define("Multijet_B", "RawPuppiMET_polar.Dot(Multijet_leadPolar) / (2.0 * Multijet_ptAvg * Multijet_leadPolar.R())")
                    .Define("Multijet_Brecoil", "RawPuppiMET_polar.Dot(Multijet_recoilPolar) / (2.0 * Multijet_ptAvg * Multijet_recoilPolar.R())")
                    .Define("Multijet_mpfResponse", "(1.0 + Multijet_B) / (1.0 - Multijet_B)")
                    .Define("Multijet_mpfLeadResponse", "1.0 - RawPuppiMET_polar.Dot(Multijet_leadPolar) / (Multijet_leadPolar.R() * Multijet_leadPolar.R())")
                    .Define("Multijet_mpfAvgResponse", "1.0 + RawPuppiMET_polar.Dot(Multijet_avgPolar) / (Multijet_avgPolar.R() * Multijet_avgPolar.R())")
                    .Define("Multijet_mpfAvpResponse", "1.0 + RawPuppiMET_polar.Dot(Multijet_bisectorPolar) / (Multijet_bisectorPolar.R() * Multijet_bisectorPolar.R())")
                    .Define("Multijet_mpfRecoilResponse", "1.0 - RawPuppiMET_polar.Dot(Multijet_recoilPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")

                    )

        return rdf_multijet
        