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
                db_rdf.Histo2D((f"DB_{system}_PtAvgVsA", "DB_"+ str(system) + "_PtVsAsymmetry;p_{T, ave} (GeV);asymmetry;N_{events}",
                                self.bins["pt"]["n"], self.bins["pt"]["bins"], self.bins["asymmetry"]["n"], self.bins["asymmetry"]["bins"]),
                                "Dijet_avgPt", "Dijet_A", "weight"),
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
                # Vs Run
                mpf_rdf.Profile1D((f"MPF_{system}_RunVsResponse", "MPF_"+ str(system) + "_RunVsResponse;Run;response;N_{events}",
                                self.bins["runs"]["n"], self.bins["runs"]["bins"]),
                                "run", "Dijet_mpfResponseCorrected", "weight"),
            ])

        return self
    
    
    def __sample_cut(self, rdf : RNode) -> RNode:
        min_pt = 15
        tag_eta = 1.3
        asymmetry_alpha = 0.7
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
                    .Define("Dijet_activityEta", "Jet_eta[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityPhi", "Jet_phi[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityMass", "Jet_mass[Jet_order != tag_idx && Jet_order != probe_idx]")
                    .Define("Dijet_activityVector", "sum_as_four_vectors(Dijet_activityPt, Dijet_activityEta, Dijet_activityPhi, Dijet_activityMass)")
                    .Define("Dijet_activityPolar", "ROOT::Math::Polar2DVectorF(Dijet_activityVector.Pt(), Dijet_activityVector.Phi())")
                    # Vectors for the tag and probe jets
                    .Define("Dijet_tagPolar", "ROOT::Math::Polar2DVectorF(Jet_pt[tag_idx], Jet_phi[tag_idx])")
                    .Define("Dijet_probePolar", "ROOT::Math::Polar2DVectorF(Jet_pt[probe_idx], Jet_phi[probe_idx])")
                    .Define("Dijet_bisectorPolar", "ROOT::Math::Polar2DVectorF(1.0, (Dijet_tagPolar - Dijet_probePolar).Phi())")
                    .Define("Dijet_metPolar", "ROOT::Math::Polar2DVectorF(RawPuppiMET_pt, RawPuppiMET_phi)")
                    .Define("Dijet_unclusteredPolar", "Dijet_metPolar + Dijet_activityPolar + Dijet_tagPolar + Dijet_probePolar")
                    .Define("deltaEta_dijet", "abs(Jet_eta[tag_idx] - Jet_eta[probe_idx])")
                    .Define("deltaR_dijet", "sqrt(deltaPhi_dijet*deltaPhi_dijet + deltaEta_dijet*deltaEta_dijet)")
                    .Define("Dijet_avgPt", "(Jet_pt[tag_idx] + Jet_pt[probe_idx]) / 2.0")
                    # Filter the two jets
                    .Filter(f"(abs(Jet_pt[tag_idx] - Jet_pt[probe_idx]) / (2.0 * Dijet_avgPt)) < {asymmetry_alpha}")
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
                    .Define("Dijet_mpfResponseCorrected", "Dijet_mpfResponse - (Dijet_tagPolar.Dot(Dijet_activityPolar + Dijet_unclusteredPolar)) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())")
                    .Define("Dijet_angleCorrection", "-ROOT::VecOps::DeltaPhi(Dijet_tagPolar.Phi(), Dijet_probePolar.Phi())")
                    # Responses as defined in the previous Dijet code
                    .Define("Dijet_mpfTagResponse", "1.0 + (Dijet_metPolar.Dot(Dijet_tagPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R())) ")
                    .Define("Dijet_dbTagResponse", "1.0 + ((-Dijet_tagPolar - Dijet_probePolar).Dot(Dijet_tagPolar) / (Dijet_tagPolar.R() * Dijet_tagPolar.R()))")
                    .Define("Dijet_mpfProbeResponse", "1.0 + (Dijet_metPolar.Dot(-Dijet_probePolar) / (Dijet_probePolar.R() * Dijet_probePolar.R()))")
                    .Define("Dijet_dbProbeResponse", "1.0 + ((-Dijet_tagPolar - Dijet_probePolar).Dot(-Dijet_probePolar) / (Dijet_probePolar.R() * Dijet_probePolar.R()))")
                    .Define("Dijet_mpfAvgResponse", "1.0 + (Dijet_metPolar.Dot(Dijet_bisectorPolar) / (Dijet_bisectorPolar.R() * Dijet_avgPt))")
                    .Define("Dijet_dbAvgResponse", "1.0 + ((-Dijet_tagPolar - Dijet_probePolar).Dot(Dijet_bisectorPolar) / (Dijet_bisectorPolar.R() * Dijet_avgPt))")
                    )
        return rdf_dijet
    
        