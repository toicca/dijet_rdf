import ROOT
from typing import List, Dict
import numpy as np
from RDFAnalyzer import RDFAnalyzer
    
RDataFrame = ROOT.RDataFrame
RunGraphs = ROOT.RDF.RunGraphs
RNode = ROOT.RDF.RNode
    
class DijetAnalyzer(RDFAnalyzer):
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
                selection_only : bool = True,
                header_dir : str = "src"
                ):
        super().__init__(filelist, trigger_list, trigger_details, json_file, nFiles, JEC, nThreads, progress_bar, isMC=isMC, local=local, system="dijet", run_raw=run_raw, selection_only=selection_only, header_dir=header_dir)
        
    def Flag_cut(self, rdf: RNode) -> RNode:
        return super().Flag_cut(rdf)
 
    def do_sample_control(self) -> "DijetAnalyzer":
        print(f"Creating control histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            control_rdf = self.__sample_cut(self.Flag_cut(rdf))
            
            self.histograms[trigger].extend([
                control_rdf.Histo1D((f"Control_{self.system}_PtAvg", "Control_"+ str(self.system) + "_PtAvg;p_{T, ave} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_avgPt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_PtProbe", "Control_"+ str(self.system) + "_PtProbe;p_{T, probe} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_probePt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_PtTag", "Control_"+ str(self.system) + "_PtTag;p_{T, tag} (GeV);N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_tagPt", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_EtaProbe", "Control_"+ str(self.system) + "_EtaProbe;#eta_{probe};N_{events}",
                                    self.bins["eta"]["n"], self.bins["eta"]["bins"]),
                                    "Dijet_probeEta", "weight"),
                control_rdf.Histo1D((f"Control_{self.system}_deltaPhi", "Control_"+ str(self.system) + "_deltaPhi;#Delta#phi_{dijet};N_{events}",
                                    self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                    "deltaPhi_dijet", "weight"),
                control_rdf.Profile1D((f"Control_{self.system}_PtAvgVsTransversePt", "Control_"+ str(self.system) + "_PtAvgVsTransversePt;p_{T, avg} (GeV);transverse p_{T};N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_avgPt", "Dijet_transversePt", "weight"),
                # Direct response as a function of pT
                control_rdf.Profile1D((f"Control_{self.system}_angCorrectionVsPt", "Control_"+ str(self.system) + "_angCorrectionVsPt;p_{T, avg} (GeV);angle correction;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                    "Dijet_avgPt", "Dijet_angleCorrection", "weight"),
            ])
            
        return self
 

    def do_DB(self) -> "DijetAnalyzer":
        system = self.system
        print(f"Creating DB histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.__sample_cut(self.Flag_cut(rdf))    
            
            self.histograms[trigger].extend([
                # 3D Distributions from which information can be projected out
                db_rdf.Histo3D((f"DB_{system}_PtTagVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, tag} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_tagPt", "Dijet_probeEta", "Dijet_dbResponseCorrected", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtTagVsEtaVsTagResponse", "MPF_"+ str(system) + "_PtVsEtaVsTagResponse;p_{T, tag} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_tagPt", "Dijet_probeEta", "Dijet_dbTagResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtProbeVsEtaVsProbeResponse", "MPF_"+ str(system) + "_PtVsEtaVsProbeResponse;p_{T, probe} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_probePt", "Dijet_probeEta", "Dijet_dbProbeResponse", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtAvgVsEtaVsAvgResponse", "MPF_"+ str(system) + "_PtVsEtaVsAvgResponse;p_{T, avg} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_avgPt", "Dijet_probeEta", "Dijet_dbAvgResponse", "weight"),
                # 3D asymmetry distribution for JER 
                db_rdf.Histo3D((f"DB_{system}_PtAvgVsEtaVsA", "DB_"+ str(system) + "_PtVsEtaVsAsymmetry;p_{T, avg} (GeV);#eta_{probe};asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Dijet_avgPt", "Dijet_probeEta", "Dijet_A", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtProbeVsEtaVsA", "DB_"+ str(system) + "_PtVsEtaVsAsymmetry;p_{T, probe} (GeV);#eta_{probe};asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Dijet_probePt", "Dijet_probeEta", "Dijet_A", "weight"),
                db_rdf.Histo3D((f"DB_{system}_PtTagVsEtaVsA", "DB_"+ str(system) + "_PtVsEtaVsAsymmetry;p_{T, tag} (GeV);#eta_{probe};asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Dijet_tagPt", "Dijet_probeEta", "Dijet_A", "weight"),
                db_rdf.Histo2D((f"DB_{system}_PtAvgVsA", "DB_"+ str(system) + "_PtVsAsymmetry;p_{T, ave} (GeV);asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Dijet_avgPt", "Dijet_A", "weight"),
                db_rdf.Histo2D((f"DB_{system}_PtProbeVsA", "DB_"+ str(system) + "_PtVsAsymmetry;p_{T, probe} (GeV);asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Dijet_probePt", "Dijet_A", "weight"),
                db_rdf.Histo2D((f"DB_{system}_PtTagVsA", "DB_"+ str(system) + "_PtVsAsymmetry;p_{T, tag} (GeV);asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Dijet_tagPt", "Dijet_A", "weight"),
                # 2D Asymmetry histogram for veto maps
                db_rdf.Profile2D((f"DB_{system}_EtaprobeVsPhiprobeVsAsymmetry", "DB_"+ str(system) + "_EtaprobeVsPhiprobeVsAsymmetry;#eta_{probe};#phi_{probe};asymmetry;N_{events}",
                                  self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["phi"]["n"], self.bins["phi"]["bins"]),
                                  "Dijet_probeEta", "Dijet_probePhi", "Dijet_A", "weight"),
                # 1D Distributions for quick comparison
                db_rdf.Profile1D((f"DB_{system}_PtTagVsResponse", "DB_"+ str(system) + "_PtTagVsResponse;p_{T, tag} (GeV);response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_dbResponse", "weight"),
                db_rdf.Profile1D((f"DB_{system}_PtTagVsResponseCorrected", "DB_"+ str(system) + "_PtTagVsResponse;p_{T, tag} (GeV);response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_dbResponseCorrected", "weight"),
            ])

            if not self.isMC:
                self.histograms[trigger].extend([
                    # Vs Run
                    db_rdf.Profile1D((f"DB_{system}_RunVsResponse", "DB_"+ str(system) + "_RunVsResponse;Run;response;N_{events}",
                                    self.bins["runs"]["n"], self.bins["runs"]["bins"]),
                                    "run", "Dijet_dbResponseCorrected", "weight"),
                ])
            barrel_rdf = db_rdf.Filter("abs(Dijet_probeEta) < 1.3")

            self.histograms[trigger].extend([
                barrel_rdf.Histo2D((f"DB_{system}_PtAvgVsA_Barrel", "DB_"+ str(system) + "_PtAvgVsA_Barrel;p_{T, avg} (GeV);asymmetry;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                    "Dijet_avgPt", "Dijet_A", "weight"),
                barrel_rdf.Histo2D((f"DB_{system}_PtTagVsA_Barrel", "DB_"+ str(system) + "_PtTagVsA_Barrel;p_{T, tag} (GeV);response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                    "Dijet_tagPt", "Dijet_A", "weight"),
                barrel_rdf.Histo2D((f"DB_{system}_PtTagVsA_Barrel", "DB_"+ str(system) + "_PtTagVsA_Barrel;p_{T, tag} (GeV);response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                    "Dijet_tagPt", "Dijet_A", "weight"),
            ])

        return self

    
    def do_MPF(self) -> "DijetAnalyzer":
        system = self.system
        print(f"Creating MPF histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            mpf_rdf = self.__sample_cut(self.Flag_cut(rdf))
            
            self.histograms[trigger].extend([
                # 3D Distributions
                mpf_rdf.Histo3D((f"MPF_{system}_PtTagVsEtaVsResponse", "MPF_"+ str(system) + "_PtVsEtaVsResponse;p_{T, tag} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_tagPt", "Dijet_probeEta", "Dijet_mpfResponseCorrected", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtTagVsEtaVsTagResponse", "MPF_"+ str(system) + "_PtVsEtaVsTagResponse;p_{T, tag} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_tagPt", "Dijet_probeEta", "Dijet_mpfTagResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtProbeVsEtaVsProbeResponse", "MPF_"+ str(system) + "_PtVsEtaVsProbeResponse;p_{T, probe} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_probePt", "Dijet_probeEta", "Dijet_mpfProbeResponse", "weight"),
                mpf_rdf.Histo3D((f"MPF_{system}_PtAvgVsEtaVsAvgResponse", "MPF_"+ str(system) + "_PtVsEtaVsAvgResponse;p_{T, avg} (GeV);#eta_{probe};response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_avgPt", "Dijet_probeEta", "Dijet_mpfAvgResponse", "weight"),
                # 1D Distributions
                mpf_rdf.Profile1D((f"MPF_{system}_PtTagVsResponse", "MPF_"+ str(system) + "_PtTagVsResponse;p_{T, tag} (GeV);response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_mpfResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtTagVsResponseCorrected", "MPF_"+ str(system) + "_PtTagVsResponse;p_{T, tag} (GeV);response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_mpfResponseCorrected", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtTagVsActivityResponse", "MPF_"+ str(system) + "_PtTagVsActivityResponse;p_{T, tag} (GeV);response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_mpfActivityResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtTagVsUnclusteredResponse", "MPF_"+ str(system) + "_PtTagVsUnclusteredResponse;p_{T, tag} (GeV);response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_mpfUnclusteredResponse", "weight"),
                mpf_rdf.Profile1D((f"MPF_{system}_PtTagVsRawMETResponse", "MPF_"+ str(system) + "_PtTagVsRawMETResponse;p_{T, tag} (GeV);response;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_mpfRawMETResponse", "weight"),                                
            ])
            if not self.isMC:
                self.histograms[trigger].extend([
                # Vs Run
                mpf_rdf.Profile1D((f"MPF_{system}_RunVsResponse", "MPF_"+ str(system) + "_RunVsResponse;Run;response;N_{events}",
                                self.bins["runs"]["n"], self.bins["runs"]["bins"]),
                                "run", "Dijet_mpfResponse", "weight"),
                ])


        return self

    def do_Activity(self) -> "DijetAnalyzer":
        system = self.system
        print(f"Creating resolution histograms for system: {self.system}")
        for trigger, rdf in self.trigger_rdfs.items():
            res_rdf = self.__sample_cut(self.Flag_cut(rdf)).Filter("Jet_pt[0] > 400").Filter("Dijet_probeEta < 1.3")
            
            self.histograms[trigger].extend([
                res_rdf.Histo1D((f"Activity_{system}_AId", "Activity_"+ str(system) + "_Id;JetId;N_{events}",
                                7, 0, 7),
                                "Dijet_activityId", "weight"),
                res_rdf.Histo1D((f"Activity_{system}_JetId", "Activity_"+ str(system) + "_JetId;JetId;N_{events}",
                                7, 0, 7),
                                "Jet_jetId", "weight"),
                res_rdf.Histo2D((f"Activity_{system}_APtVsMPF", "Activity_"+ str(system) + "_Pt;p_{T,activity} (GeV);MPF",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_activityPtLen", "Dijet_mpfResponse", "weight"),
                res_rdf.Histo2D((f"Activity_{system}_IdSumVsMPF", "Activity_"+ str(system) + "_IdSum;IdSum;MPF",
                                 20, 0, 20, self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_activityIdSum", "Dijet_mpfResponse", "weight"),
                res_rdf.Profile1D((f"Activity_{system}_ProfileAPtVsMPF", "Activity_"+ str(system) + "_Pt;p_{T,activity} (GeV);MPF",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_activityPtLen", "Dijet_mpfResponse", "weight"),
                res_rdf.Profile1D((f"Activity_{system}_ProfileIdSumVsMPF", "Activity_"+ str(system) + "_IdSum;Average JetId;MPF",
                                7, 0, 7),
                                "Dijet_activityIdSum", "Dijet_mpfResponse", "weight"),
                res_rdf.Histo1D((f"Activity_{system}_MPF", "Activity_"+ str(system) + "_MPF;MPF;N_{events}",
                                self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_mpfResponse", "weight"),
                res_rdf.Histo1D((f"Activity_{system}_MPFCorrected", "Activity_"+ str(system) + "_MPF;MPF;N_{events}",
                                self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_mpfResponseCorrected", "weight"),
                res_rdf.Histo2D((f"Activity_{system}_PtAvgVsPtActivity", "Activity_"+ str(system) + "_PtAvgVsPtActivity;p_{T, avg} (GeV);p_{T, activity} (GeV)",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_avgPt", "Dijet_activityPtLen", "weight"),
                res_rdf.Profile1D((f"Activity_{system}_ProfilePtAvgVsPtActivity", "Activity_"+ str(system) + "_PtAvgVsPtActivity;p_{T, avg} (GeV);p_{T, activity} (GeV)",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_avgPt", "Dijet_activityPtLen", "weight"),
                res_rdf.Profile1D((f"Activity_{system}_PtActivityVsActivityResponse", "Activity_"+ str(system) + "_PtActivityVsActivityResponse;p_{T, activity} (GeV);MPF",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_activityPtLen", "Dijet_activityResponse", "weight"),
                res_rdf.Profile1D((f"Activity_{system}_PtTagVsActivityResponse", "Activity_"+ str(system) + "_PtTagVsActivityResponse;p_{T, tag} (GeV);MPF",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_tagPt", "Dijet_activityResponse", "weight"),
            ])


            cut_id_rdf = res_rdf.Filter("Dijet_activityIdSum > 3")

            self.histograms[trigger].extend([
                cut_id_rdf.Histo1D((f"Activity_{system}_AId_Cut", "Activity_"+ str(system) + "_Id;JetId;N_{events}",
                                7, 0, 7),
                                "Dijet_activityId", "weight"),
                cut_id_rdf.Histo2D((f"Activity_{system}_APtVsMPF_Cut", "Activity_"+ str(system) + "_Pt;p_{T,activity} (GeV);MPF",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_activityPtLen", "Dijet_mpfResponse", "weight"),
                cut_id_rdf.Histo2D((f"Activity_{system}_IdSumVsMPF_Cut", "Activity_"+ str(system) + "_IdSum;IdSum;MPF",
                                 20, 0, 20, self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_activityIdSum", "Dijet_mpfResponse", "weight"),
                cut_id_rdf.Profile1D((f"Activity_{system}_ProfileAPtVsMPF_Cut", "Activity_"+ str(system) + "_Pt;p_{T,activity} (GeV);MPF",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"]),
                                "Dijet_activityPtLen", "Dijet_mpfResponse", "weight"),
                cut_id_rdf.Profile1D((f"Activity_{system}_ProfileIdSumVsMPF_Cut", "Activity_"+ str(system) + "_IdSum;Average JetId;MPF",
                                7, 0, 7),
                                "Dijet_activityIdSum", "Dijet_mpfResponse", "weight"),
                cut_id_rdf.Histo1D((f"Activity_{system}_MPF_Cut", "Activity_"+ str(system) + "_MPF;MPF;N_{events}",
                                self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_mpfResponse", "weight"),
                cut_id_rdf.Histo1D((f"Activity_{system}_MPFCorrected_Cut", "Activity_"+ str(system) + "_MPF;MPF;N_{events}",
                                self.bins["response"]["n"], self.bins["response"]["bins"]),
                                "Dijet_mpfResponseCorrected", "weight"),
            ])
                                    

        return self
    
    
    def __sample_cut(self, rdf : RNode) -> RNode:
        min_pt = 0 #15
        tag_eta = 1.3
        asymmetry_alpha = 1.0
        thirdJet_alpha = 0.2
        delta_phi = 2.7

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

        rdf_dijet = (rdf.Filter("nJet >= 2")
                    # Choose tag and probe jets, use gRandom-Rndm(), should be made deterministic
                    .Define("tag_idx", "int(round(gRandom->Rndm()))")
                    .Define("probe_idx", "tag_idx == 0")
                    .Filter(f"Jet_pt[tag_idx] > {min_pt} && Jet_pt[probe_idx] > {min_pt}")
                    .Filter("Jet_jetId[tag_idx] >= 4 && Jet_jetId[probe_idx] >= 4")
                    .Define("deltaPhi_dijet", "abs(Jet_phi[tag_idx] - Jet_phi[probe_idx])")
                    .Redefine("deltaPhi_dijet", "deltaPhi_dijet > TMath::Pi() ? TMath::TwoPi() - deltaPhi_dijet : deltaPhi_dijet")
                    .Filter(f"abs(Jet_eta[tag_idx]) < {tag_eta} && deltaPhi_dijet > {delta_phi}")
                    .Define("Dijet_activityPt", "Jet_pt[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityPtSum", "ROOT::VecOps::Sum(Dijet_activityPt)")
                    .Filter("(Jet_pt[tag_idx] + Jet_pt[probe_idx]) > Dijet_activityPtSum * 2")
                    .Define("Dijet_activityEta", "Jet_eta[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityPhi", "Jet_phi[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityMass", "Jet_mass[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityId", "Jet_jetId[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityIdSum", "float(ROOT::VecOps::Sum(Dijet_activityId)) / float(Dijet_activityId.size())")
                    .Define("Dijet_activityVector", "sum_as_four_vectors(Dijet_activityPt, Dijet_activityEta, Dijet_activityPhi, Dijet_activityMass)")
                    .Define("Dijet_activityPolar", "ROOT::Math::Polar2DVectorF(Dijet_activityVector.Pt(), Dijet_activityVector.Phi())")
                    .Define("Dijet_activityPtLen", "Dijet_activityVector.Pt()")
                    .Define("Dijet_eventVector", f"sum_as_four_vectors(Jet_pt[Jet_pt >= {min_pt}], Jet_eta[Jet_pt >= {min_pt}], Jet_phi[Jet_pt >= {min_pt}], Jet_mass[Jet_pt >= {min_pt}])")
                    .Define("Dijet_eventPolar", "ROOT::Math::Polar2DVectorF(Dijet_eventVector.Pt(), Dijet_eventVector.Phi())")
                    .Define("Dijet_rawVector", f"sum_as_four_vectors(Jet_pt[Jet_pt >= {min_pt}] * (1.0 - Jet_rawFactor[Jet_pt >= {min_pt}]), Jet_eta[Jet_pt >= {min_pt}], Jet_phi[Jet_pt >= {min_pt}], Jet_mass[Jet_pt >= {min_pt}] * (1.0 - Jet_rawFactor[Jet_pt >= {min_pt}]))")
                    .Define("Dijet_rawPolar", "ROOT::Math::Polar2DVectorF(Dijet_rawVector.Pt(), Dijet_rawVector.Phi())")
                    .Define("Dijet_t1METVector", "Dijet_rawVector - Dijet_eventVector")
                    .Define("Dijet_t1METPolar", "ROOT::Math::Polar2DVectorF(Dijet_t1METVector.Pt(), Dijet_t1METVector.Phi())")
                    # Vectors for the tag and probe jets
                    .Define("Dijet_tagPolar", "ROOT::Math::Polar2DVectorF(Jet_pt[tag_idx], Jet_phi[tag_idx])")
                    .Define("Dijet_probePolar", "ROOT::Math::Polar2DVectorF(Jet_pt[probe_idx], Jet_phi[probe_idx])")
                    .Define("Dijet_bisectorPolar", "ROOT::Math::Polar2DVectorF(1.0, (Dijet_tagPolar - Dijet_probePolar).Phi())")
                    .Define("Dijet_metPolar", "ROOT::Math::Polar2DVectorF(Dijet_t1METPolar.R(), Dijet_t1METPolar.Phi())") # Define through T1MET
                    .Define("Dijet_unclusteredPolar", "Dijet_metPolar + Dijet_activityPolar + Dijet_tagPolar + Dijet_probePolar")
                    .Define("deltaEta_dijet", "abs(Jet_eta[tag_idx] - Jet_eta[probe_idx])")
                    .Define("deltaR_dijet", "sqrt(deltaPhi_dijet*deltaPhi_dijet + deltaEta_dijet*deltaEta_dijet)")
                    .Define("Dijet_avgPt", "(Jet_pt[tag_idx] + Jet_pt[probe_idx]) / 2.0")
                    .Define("Dijet_transversePt", "(Dijet_tagPolar + Dijet_probePolar + Dijet_metPolar + Dijet_activityPolar).R()")
                    .Define("Dijet_rawPtSum", "(RawPuppiMET_polar + Dijet_rawPolar).R()")
                    # Filter the two jets
                    .Filter(f"(abs(Jet_pt[tag_idx] - Jet_pt[probe_idx]) / (2.0 * Dijet_avgPt)) < {asymmetry_alpha}")
                    .Filter(f"nJet > 2 ? Jet_pt[2] / Dijet_avgPt < {thirdJet_alpha} : true")
                    .Define("Dijet_tagPt", "Jet_pt[tag_idx]")
                    .Define("Dijet_probePt", "Jet_pt[probe_idx]")
                    .Define("Dijet_tagEta", "Jet_eta[tag_idx]") # Scuffed
                    .Define("Dijet_probeEta", "Jet_eta[probe_idx]")
                    .Define("Dijet_tagPhi", "Jet_phi[tag_idx]")
                    .Define("Dijet_probePhi", "Jet_phi[probe_idx]")
                    # Responses
                    .Define("Dijet_A", "(Jet_pt[probe_idx] - Jet_pt[tag_idx]) / (2.0 * Dijet_avgPt)")
                    .Define("Dijet_dbResponse", "Dijet_probePolar.R() / Dijet_tagPolar.R()")
                    .Define("Dijet_dbResponseCorrected", "-Dijet_probePolar.Dot(Dijet_tagPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())")
                    .Define("Dijet_mpfResponse", "1 + Dijet_metPolar.Dot(Dijet_tagPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())")
                    .Define("Dijet_mpfResponseCorrected", "Dijet_mpfResponse - (Dijet_tagPolar.Dot(Dijet_activityPolar)) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())")
                    .Define("Dijet_mpfActivityResponse", "1 - Dijet_tagPolar.Dot(Dijet_activityPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())")
                    .Define("Dijet_mpfUnclusteredResponse", "1 - Dijet_tagPolar.Dot(Dijet_unclusteredPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())")
                    .Define("Dijet_mpfRawMETResponse", "1 + Dijet_tagPolar.Dot(RawPuppiMET_polar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())")
                    .Define("Dijet_angleCorrection", "-ROOT::VecOps::DeltaPhi(Dijet_tagPolar.Phi(), Dijet_probePolar.Phi())")
                    # Responses as defined in the previous Dijet code
                    .Define("Dijet_mpfTagResponse", "1.0 + (Dijet_metPolar.Dot(Dijet_tagPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())) ")
                    .Define("Dijet_dbTagResponse", "1.0 + ((-Dijet_tagPolar - Dijet_probePolar).Dot(Dijet_tagPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R()))")
                    .Define("Dijet_mpfProbeResponse", "1.0 + (Dijet_metPolar.Dot(-Dijet_probePolar) / (Dijet_probePolar.R() * Dijet_probePolar.R()))")
                    .Define("Dijet_dbProbeResponse", "1.0 + ((-Dijet_tagPolar - Dijet_probePolar).Dot(-Dijet_probePolar) / (Dijet_probePolar.R() * Dijet_probePolar.R()))")
                    .Define("Dijet_mpfAvgResponse", "1.0 + (Dijet_metPolar.Dot(Dijet_bisectorPolar) / (Dijet_bisectorPolar.R() * Dijet_avgPt))")
                    .Define("Dijet_dbAvgResponse", "1.0 + ((-Dijet_tagPolar - Dijet_probePolar).Dot(Dijet_bisectorPolar) / (Dijet_bisectorPolar.R() * Dijet_avgPt))")
                    .Define("Dijet_activityResponse", "1.0 - (Dijet_tagPolar.Dot(Dijet_activityPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R()))")
                    )
        return rdf_dijet
    
        