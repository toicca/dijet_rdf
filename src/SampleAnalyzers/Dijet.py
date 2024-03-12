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
                ):
        super().__init__(filelist, trigger_list, json_file, nFiles, JEC, nThreads, progress_bar)
        
    def Flag_cut(self, rdf: RNode) -> RNode:
        return super().Flag_cut(rdf)
 

    def do_DB(self) -> "DijetAnalyzer":
        system = "dijet"
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.Flag_cut(self.__sample_cut(rdf))
            pT_binLabels = ["average_Pt_dijet", "Jet_pt_tag", "Jet_pt_probe"]

            db_rdf = (db_rdf.Define(f"response_DB_{system}", f"(1.0 + Asymmetry_{system}) / (1.0 - Asymmetry_{system})")
                    )
            
            print("Creating DB histograms for trigger", trigger)
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"DB_{system}_Response", "DB_"+ str(system) + "_response;response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), f"response_DB_{system}", "weight"),
                db_rdf.Histo2D((f"DB_{system}_EtaVsResponse", "DB_"+ str(system) + "_EtaVsResponse;|#eta|;response", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), f"Jet_eta_tag_{system}", f"response_DB_{system}", "weight")
            ])
            for label in pT_binLabels:
                self.histograms[trigger].append(
                    db_rdf.Histo2D((f"DB_{system}_PtVsResponse_{label}", "DB_"+str(system)+"_PtVsResponse_"+label+";p_{T, "+label+"} (GeV);response;N_{events}",
                                    self.bins["pt"]["n"], self.bins["pt"]["bins"],
                                    self.bins["response"]["n"], self.bins["response"]["bins"]),
                                    label, f"response_DB_{system}", "weight"),
                )

        return self
    
    def do_MPF(self) -> "DijetAnalyzer":
        system = "dijet"
        for trigger, rdf in self.trigger_rdfs.items():
            db_rdf = self.Flag_cut(self.__sample_cut(rdf))
            pT_binLabels = ["average_Pt_dijet", "Jet_pt_tag", "Jet_pt_probe"]

            db_rdf = (db_rdf.Define(f"response_MPF_{system}", f"(1.0 + B_{system}) / (1.0 - B_{system})")
                    )
            
            print("Creating MPF histograms for trigger", trigger)
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"MPF_{system}_Response", "MPF_"+ str(system) + "_response;response;N_{events}", self.bins["response"]["n"], self.bins["response"]["bins"]), 
                               f"response_MPF_{system}", "weight"),
                db_rdf.Histo2D((f"MPF_{system}_EtaVsResponse", "DB_"+ str(system) + "_EtaVsResponse;|#eta|;response", self.bins["eta"]["n"], self.bins["eta"]["bins"], self.bins["response"]["n"], self.bins["response"]["bins"]), 
                               f"Jet_eta_tag_{system}", f"response_MPF_{system}", "weight")
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
        print("Applying sample cuts")
        min_pt = 15
        tag_eta = 1.3
        asymmetry_alpha = 0.7
        delta_phi = 2.7
        rdf_dijet = (rdf.Filter("nJet >= 2")
                    # Choose tag and probe jets, this isn't necessarily random so check on it
                    .Define("tag_idx", "int(run % 2 == 0)")
                    .Define("probe_idx", "tag_idx == 0")
                    .Define("deltaPhi_dijet", "abs(Jet_phi[tag_idx] - Jet_phi[probe_idx])")
                    .Define("deltaEta_dijet", "abs(Jet_eta[tag_idx] - Jet_eta[probe_idx])")
                    .Define("deltaR_dijet", "sqrt(deltaPhi_dijet*deltaPhi_dijet + deltaEta_dijet*deltaEta_dijet)")
                    .Define("average_Pt_dijet", "(Jet_pt[tag_idx] + Jet_pt[probe_idx]) / 2.0")
                    .Define("Jet_pt_tag", "Jet_pt[tag_idx]")
                    .Define("Jet_pt_probe", "Jet_pt[probe_idx]")
                    .Define("Jet_eta_tag_dijet", "Jet_eta[tag_idx]") # Scuffed
                    # Filter the two jets
                    .Filter(f"abs(Jet_eta[tag_idx]) < {tag_eta} && deltaPhi_dijet > {delta_phi}")
                    .Filter(f"Jet_pt[tag_idx] > {min_pt} && Jet_pt[probe_idx] > {min_pt}")
                    .Filter("Jet_jetId[tag_idx] >= 4 && Jet_jetId[probe_idx] >= 4")
                    .Filter(f"(abs(Jet_pt[tag_idx] - Jet_pt[probe_idx]) / (2.0 * average_Pt_dijet)) < {asymmetry_alpha}")
                    # Asymmetry of the system
                    .Define("Asymmetry_dijet", "(Jet_pt[tag_idx] - Jet_pt[probe_idx]) / (Jet_pt[tag_idx] + Jet_pt[probe_idx])")
                    .Define("B_dijet", "RawPuppiMET_pt * cos(ROOT::VecOps::DeltaPhi(Jet_phi[tag_idx], RawPuppiMET_phi)) / (Jet_pt[tag_idx] + Jet_pt[probe_idx])")
                    # TODO: Add jet veto
                    )
        return rdf_dijet
    
        