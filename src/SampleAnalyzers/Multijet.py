import ROOT
from typing import List, Dict
import numpy as np
from RDFAnalyzer import RDFAnalyzer
    
RDataFrame = ROOT.RDataFrame
RunGraphs = ROOT.RDF.RunGraphs
RNode = ROOT.RDF.RNode
    
class MultijetAnalyzer(RDFAnalyzer):
    def __init__(self, filelist : List[str],
                trigger_list : List[str],
                trigger_details: bool = False,
                json_file : str = "",
                nFiles : int = -1,
                JEC : Dict = {},
                nThreads : int = 1,
                progress_bar : bool = False,
                isMC : bool = False,
                local : bool = False,
                run_raw : bool = False,
                selection_only : bool = False,
                header_dir : str = "src"
                ):
        super().__init__(filelist, trigger_list, trigger_details, json_file, nFiles, JEC, nThreads, progress_bar, isMC=isMC, local=local, system="multijet", run_raw=run_raw, selection_only=selection_only, header_dir=header_dir)
 
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
                                    
                # Unclustered energy and activity
                control_rdf.Histo1D((f"Control_{self.system}_UnclusteredPt", "Control_"+ str(self.system) + "_UnclusteredPt;p_{T, unclustered} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_unclusteredPt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_UnclusteredPhi", "Control_"+ str(self.system) + "_UnclusteredPhi;#phi_{unclustered};N_{events}",
                                    self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                    "Multijet_unclusteredPhi", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_ActivityPt", "Control_"+ str(self.system) + "_ActivityPt;p_{T, activity} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_activityPt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_ActivityPhi", "Control_"+ str(self.system) + "_ActivityPhi;#phi_{activity};N_{events}",
                                    self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                    "Multijet_activityPhi", "weight"),
                control_rdf.Profile1D((f"Control_{self.system}_PtRecoilVsActivityPt", "Control_"+ str(self.system) + "_PtRecoilVsActivityPt;p_{T, recoil} (GeV);p_{T, activity} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_recoilPt", "Multijet_activityPt", "weight"),
                control_rdf.Profile1D((f"Control_{self.system}_PtRecoilVsUnclusteredPt", "Control_"+ str(self.system) + "_PtRecoilVsUnclusteredPt;p_{T, recoil} (GeV);p_{T, unclustered} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_recoilPt", "Multijet_unclusteredPt", "weight"),
                # T1MET
                control_rdf.Profile1D((f"Control_{self.system}_PtRecoilVst1METPt", "Control_"+ str(self.system) + "_PtRecoilVst1METPt;p_{T, recoil} (GeV);p_{T, t1MET} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Multijet_recoilPt", "Multijet_t1METPt", "weight"),
        ])
            
        return self

    def do_DB(self) -> "MultijetAnalyzer":
        system = "multijet"
        print(f"Creating DB histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))    
            
            self.histograms[trigger].extend([
                # 3D Distributions
                db_rdf.Histo3D((f"DB_{system}_PtRecoilVsEtaVsResponse", "DB_"+ str(system) + "_PtVsEtaVsResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_dbResponseCorrected", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtRecoilVsEtaVsRecoilResponse", "DB_"+ str(system) + "_PtVsEtaVsRecoilResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_dbRecoilResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtLeadVsEtaVsLeadResponse", "DB_"+ str(system) + "_PtVsEtaVsLeadResponse;p_{T, lead} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_dbLeadResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtAvpVsEtaVsAvgResponse", "DB_"+ str(system) + "_PtVsEtaVsAvgResponse;p_{T, avp} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_ptAvp", "Multijet_recoilEta", "Multijet_dbAvgResponse", "weight"),
                # Average distributions for pT in barrel
                db_rdf.Profile1D((f"DB_{system}_PtRecoilVsResponseCorrected", "DB_"+ str(system) + "_PtRecoilVsResponse;p_{T, lead} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_dbResponseCorrected", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtRecoilVsResponse", "DB_"+ str(system) + "_PtRecoilVsDirectResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_dbResponse", "weight"),
            ])
            db_gt100_rdf = self.__sample_cut(self.Flag_cut(rdf)) # .Filter("Multijet_recoilPt > 100")
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
                mpf_rdf.Histo3D((f"MPF_{system}_PtRecoilVsEtaVsResponse", "MPF_"+ str(system) + "_PtRecoilVsEtaVsResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_mpfResponseCorrected", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtRecoilVsEtaVsRecoilResponse", "MPF_"+ str(system) + "_PtRecoilVsEtaVsRecoilResponse;p_{T, recoil} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_recoilPt", "Multijet_recoilEta", "Multijet_mpfRecoilResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtLeadVsEtaVsLeadResponse", "MPF_"+ str(system) + "_PtLeadVsEtaVsLeadResponse;p_{T, lead} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_leadPt", "Multijet_recoilEta", "Multijet_mpfLeadResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtAvpVsEtaVsAvgResponse", "MPF_"+ str(system) + "_PtAvpVsEtaVsAvgResponse;p_{T, avp} (GeV);#eta_{recoil};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Multijet_ptAvp", "Multijet_recoilEta", "Multijet_mpfAvgResponse", "weight"),

                # Average distributions for pT in barrel
                mpf_rdf.Profile1D((f"MPF_{system}_PtRecoilVsResponse", "MPF_"+ str(system) + "_PtRecoilVsResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_mpfResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtRecoilVsResponseCorrected", "MPF_"+ str(system) + "_PtRecoilVsResponseCorrected;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_mpfResponseCorrected", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtRecoilVsActivityResponse", "MPF_"+ str(system) + "_PtRecoilVsActivityResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_mpfActivityResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtRecoilVsUnclusteredResponse", "MPF_"+ str(system) + "_PtRecoilVsUnclusteredResponse;p_{T, recoil} (GeV);response;N_{events}",
                                        self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                        "Multijet_recoilPt", "Multijet_mpfUnclusteredResponse", "weight"),
            ])
            mpf_gt100_rdf = self.__sample_cut(self.Flag_cut(rdf)).Filter("Multijet_recoilPt > 100")
            self.histograms[trigger].extend([
                mpf_gt100_rdf.Profile1D((f"MPF_{system}_RunVsResponse", "MPF_"+ str(system) + "_RunVsResponse;Run;response;N_{events}",
                                         self.bins["runs"]["n"], self.bins["runs"]["bins"]),
                                            "run", "Multijet_mpfResponse", "weight"),
           ])

        return self
    
    def __sample_cut(self, rdf : RNode) -> RNode:
        evt_min_pt = 15
        min_pt = 30
        lead_eta = 1.3
        recoil_eta = 2.5
        asymmetry_alpha = 0.7
        delta_phi = 1.0
        activity_alpha = 0.3
        
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

        if not hasattr(ROOT, "DeltaR_to_lead"):
            ROOT.gInterpreter.Declare(""" 
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
                    .Define("Jet_pt_activity", f"Jet_pt[(Jet_order != lead_idx) && ((Jet_pt < {min_pt}) || (Jet_eta > {recoil_eta})) && (Jet_pt >= {evt_min_pt})]")
                    .Define("Jet_eta_activity", f"Jet_eta[(Jet_order != lead_idx) && ((Jet_pt < {min_pt}) || Jet_eta > {recoil_eta}) && (Jet_pt >= {evt_min_pt})]")
                    .Define("Jet_mass_activity", f"Jet_mass[(Jet_order != lead_idx) && ((Jet_pt < {min_pt}) || Jet_eta > {recoil_eta}) && (Jet_pt >= {evt_min_pt})]")
                    .Define("Jet_phi_activity", f"Jet_phi[(Jet_order != lead_idx) && ((Jet_pt < {min_pt}) || Jet_eta > {recoil_eta}) && (Jet_pt >= {evt_min_pt})]")
                    .Define("Multijet_activityVector", "sum_as_four_vectors(Jet_pt_activity, Jet_eta_activity, Jet_phi_activity, Jet_mass_activity)")
                    .Define("Multijet_recoilVector", "sum_as_four_vectors(Jet_pt_recoil, Jet_eta_recoil, Jet_phi_recoil, Jet_mass_recoil)")
                    .Define("Multijet_eventVector", f"sum_as_four_vectors(Jet_pt[Jet_pt >= {evt_min_pt}], Jet_eta[Jet_pt >= {evt_min_pt}], Jet_phi[Jet_pt >= {evt_min_pt}], Jet_mass[Jet_pt >= {evt_min_pt}])")
                    .Define("Multijet_rawVector", f"sum_as_four_vectors(Jet_pt[Jet_pt >= {evt_min_pt}] * (1.0 - Jet_rawFactor[Jet_pt >= {evt_min_pt}]), Jet_eta[Jet_pt >= {evt_min_pt}], Jet_phi[Jet_pt >= {evt_min_pt}], Jet_mass[Jet_pt >= {evt_min_pt}] * (1.0 - Jet_rawFactor[Jet_pt >= {evt_min_pt}]))")
                    .Define("Multijet_t1METVector", "Multijet_rawVector - Multijet_eventVector")
                    .Define("Multijet_t1METPolar", "ROOT::Math::Polar2DVectorF(Multijet_t1METVector.Pt(), Multijet_t1METVector.Phi())")
                    .Define("Multijet_t1METPt", "Multijet_t1METVector.Pt()")
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
                    .Define("Multijet_activityPolar", "ROOT::Math::Polar2DVectorF(Multijet_activityVector.Pt(), Multijet_activityVector.Phi())")
                    .Define("Multijet_activityPt", "Multijet_activityVector.Pt()")
                    .Define("Multijet_activityPhi", "Multijet_activityVector.Phi()")
                    .Define("Multijet_unclusteredPolar", "Multijet_recoilPolar + Multijet_activityPolar - Multijet_t1METPolar + Multijet_leadPolar")
                    .Define("Multijet_unclusteredPt", "Multijet_unclusteredPolar.R()")
                    .Define("Multijet_unclusteredPhi", "Multijet_unclusteredPolar.Phi()")
                    .Define("Multijet_avgPolar", "(Multijet_recoilPolar - Multijet_leadPolar) / 2.0")
                    .Define("Multijet_bisectorPolar", "ROOT::Math::Polar2DVectorF(1.0, (Multijet_recoilPolar / (Multijet_recoilPolar.R()) - Multijet_leadPolar / (Multijet_leadPolar.R())).Phi())")
                    .Define("Multijet_nRecoil", "Jet_pt_recoil.size()")
                    .Define("Multijet_ptAvg", "(Multijet_avgPolar.R())")
                    .Define("Multijet_ptAvp", "0.5 * (Multijet_recoilPolar.Dot(Multijet_bisectorPolar) - Multijet_leadPolar.Dot(Multijet_bisectorPolar))")
                    .Define("Multijet_deltaPhi", "ROOT::VecOps::DeltaPhi(Multijet_leadVector.Phi(), Multijet_recoilVector.Phi())")
                    .Filter(f"abs(Multijet_deltaPhi - TMath::Pi()) < 0.3")
                    .Filter(f"Multijet_activityPt < {activity_alpha} * Multijet_recoilPt")
                    .Define("Multijet_METdotRecoil", "Multijet_t1METPolar.Dot(Multijet_recoilPolar) / (Multijet_recoilPolar.R() * Multijet_t1METPolar.R())")
                    .Define("Multijet_METdotLead", "Multijet_t1METPolar.Dot(Multijet_leadPolar) / (Multijet_leadPolar.R() * Multijet_t1METPolar.R())")
                    # Responses
                    # DB
                    .Define("Multijet_A", "(Multijet_recoilPolar.R() - Multijet_leadPolar.R()) / (Multijet_recoilPolar.R() + Multijet_leadPolar.R())")
                    .Define("Multijet_dbResponseCorrected", "-Multijet_recoilPolar.Dot(Multijet_leadPolar) / (Multijet_recoilVector.Pt() * Multijet_recoilVector.Pt())") 
                    .Define("Multijet_dbResponse", "Multijet_leadVector.Pt() / Multijet_recoilVector.Pt()")
                    # MPF
                    .Define("Multijet_mpfResponse", "1.0 + Multijet_recoilPolar.Dot(Multijet_t1METPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")
                    .Define("Multijet_mpfResponseCorrected", "Multijet_mpfResponse - Multijet_recoilPolar.Dot(Multijet_activityPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")
                    .Define("Multijet_mpfActivityResponse", "1.0 - Multijet_recoilPolar.Dot(Multijet_activityPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")
                    .Define("Multijet_mpfUnclusteredResponse", "1.0 - Multijet_recoilPolar.Dot(Multijet_unclusteredPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")
                    # Responses from the previous dijet code
                    .Define("Multijet_p4m3", "-Multijet_leadPolar - Multijet_recoilPolar")
                    .Define("Multijet_dbLeadResponse", "1.0 - Multijet_p4m3.Dot(Multijet_leadPolar) / (Multijet_leadPolar.R() * Multijet_leadPolar.R())")
                    .Define("Multijet_dbAvgResponse", "1.0 + Multijet_p4m3.Dot(Multijet_bisectorPolar) / (Multijet_ptAvg)")
                    .Define("Multijet_dbRecoilResponse", "1.0 + Multijet_p4m3.Dot(Multijet_recoilPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")
                    .Define("Multijet_mpfLeadResponse", "1.0 + Multijet_t1METPolar.Dot(Multijet_leadPolar) / (Multijet_leadPolar.R() * Multijet_leadPolar.R())")
                    .Define("Multijet_mpfAvgResponse", "1.0 + Multijet_t1METPolar.Dot(Multijet_bisectorPolar) / (Multijet_ptAvg)")
                    .Define("Multijet_mpfRecoilResponse", "1.0 + Multijet_t1METPolar.Dot(Multijet_recoilPolar) / (Multijet_recoilPolar.R() * Multijet_recoilPolar.R())")
                    )

        return rdf_multijet
        
