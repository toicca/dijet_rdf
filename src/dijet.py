import ROOT
from dataclasses import dataclass
from typing import List
from bins import bins_pT, bins_eta, bins_pT_n, bins_eta_n
import numpy as np
from numba import njit
from make_JEC import compile_JEC, load_JEC, clean_JEC
    
RDataFrame = ROOT.RDataFrame
RunGraphs = ROOT.RDF.RunGraphs
    
@dataclass
class JEC_corrections:
    L1 : str
    L2Relative : str
    L2L3 : str
    
class dijet:
    def __init__(self, filelist : List[str],
                trigger_list : List[str],
                json_file : str,
                nFiles : int = -1,
                JEC : JEC_corrections = JEC_corrections("", "", ""),
                nThreads : int = 1,
                ):
        self.nThreads = nThreads
        self.trigger_list = trigger_list
        self.histograms = {"all" : []} # format : {trigger : [histograms]}
        self.trigger_rdfs = {} # format : {trigger : rdf}. self.rdf is not initialized here due to order of operations
        self.has_run = False

        ROOT.gInterpreter.Declare('#include "src/JSONRDF_code.h"')
        
        self.rdf = self.__loadRDF(filelist, nFiles = nFiles)
        ROOT.RDF.Experimental.AddProgressBar(self.rdf)

        # Initial variables
        self.rdf = (self.rdf.Define("weight", "1.0")
                    .Define("leading_idx", "0")
                    .Define("Jet_rawPt", "Jet_pt * (1.0-Jet_rawFactor)")
                    .Define("RawPuppiMET_polar", "ROOT::Math::Polar2DVectorF(RawPuppiMET_pt, RawPuppiMET_phi)")
                    )
        
        if JEC.L1 != "" or JEC.L2Relative != "" or JEC.L2L3 != "":
            # clean_JEC()
            # compile_JEC()
            load_JEC()
            self.rdf = self.__redo_JEC(JEC)

        self.rdf = self.__do_cut_golden_json(json_file)

        for trigger in trigger_list:
            # Consider changing these to a static size or Pythons array
            self.trigger_rdfs[trigger] = self.rdf.Filter(trigger)
            self.histograms[trigger] = []

        # Only after JECs, common cuts etc. have been applied, we can create the inclusive rdf
        self.trigger_rdfs["all"] = self.rdf


    def __Flag_cut(self, rdf : RDataFrame) -> RDataFrame:
        flag = """Flag_goodVertices && 
                    Flag_globalSuperTightHalo2016Filter &&
                    Flag_EcalDeadCellTriggerPrimitiveFilter &&
                    Flag_BadPFMuonFilter &&
                    Flag_BadPFMuonDzFilter && 
                    Flag_hfNoisyHitsFilter &&
                    Flag_eeBadScFilter &&
                    Flag_ecalBadCalibFilter
                    """
        rdf = (rdf.Filter(flag))
        return rdf
        
    def __loadRDF(self, filelist : List[str], treename : str = "Events", nFiles : int = -1, local : bool = False) -> RDataFrame:
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
        
        # TChain needs to be created as global, as it will otherwise go out of scope
        # With multiple analyses this might raise an issue, and the chain might need renaming
        # Another option is to use the gInterpreter, see
        # https://github.com/root-project/root/issues/10965#issue-1305124734
        global chain
        chain = ROOT.TChain(treename)
        for file in filelist[:nFiles]:
            chain.Add(file)
        
        rdf = RDataFrame(chain)
        
        return rdf
 
    def __redo_JEC(self, jec : JEC_corrections) -> RDataFrame:
        ROOT.gInterpreter.Declare('#include "src/JECRDF_code.h"')
        # Remember that this might alter the leading jets!
        ROOT.init_JEC(L1 = jec.L1, L2Relative = jec.L2Relative, L2L3 = jec.L2L3, nThreads = self.nThreads)

        rdf = (self.rdf.Define("JEC", "getJEC(rdfslot_, Jet_pt, Jet_eta, Jet_area, Rho_fixedGridRhoFastjetAll)")
               .Redefine("Jet_pt", "Jet_rawPt * JEC")
               .Redefine("Jet_rawFactor", "1.0 - 1.0 / JEC")
               .Redefine("leading_idx", "ArgMax(Jet_pt)")
        )
        self.histograms["all"].extend([rdf.Histo1D(("JEC", "JEC;JEC;N_{events}", 100, 0.0, 2.0), "JEC", "weight"),])

        return rdf
    
    def __do_cut_golden_json(self, json_file : str) -> RDataFrame:
        ROOT.init_json(json_file)
        rdf = self.rdf.Filter("isGoodLumi(run, luminosityBlock)", "JSON Filter")
        return rdf
    
    def __do_cut_veto_map(self, veto_map_file : str) -> RDataFrame:
        # TODO: Implement this
        # ROOT.init_veto_map(veto_map_file)
        rdf = self.rdf.Filter("isGoodLumi(run, luminosityBlock)", "Veto Map Filter")
        return rdf
    
    def do_smear_JER(self):
        return 0
            
    def get_histograms(self) -> dict:
        print("Returning histograms")
        return self.histograms

    def run_histograms(self):
        if not self.has_run:
            for trigger in self.trigger_list:
                print("Running histograms for trigger", trigger)
                RunGraphs(self.histograms[trigger])
            
            self.has_run = True
        else:
            print("Histograms have already been run")
            
    def do_inclusive(self):
        # Create the inclusive histograms
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = (self.__Flag_cut(rdf)).Redefine("Jet_pt", "Jet_pt[Jet_jetId >= 4]").Redefine("Jet_eta", "Jet_eta[Jet_jetId >= 4]")
            
            # Eta binned rdfs for pT distribution of jets
            eta_bins_for_pt = [(0.0, 1.3), (0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 2.5),
                                 (2.5, 3.0), (3.0, 3.5), (3.5, 4.0), (4.0, 4.5), (4.5, 5.0)]
            eta_binned_rdfs = {}
            for i, val in enumerate(eta_bins_for_pt):
                eta_binned_rdfs[i] = (selected_rdf.Redefine("Jet_pt", f"Jet_pt[abs(Jet_eta) > {val[0]} && abs(Jet_eta) < {val[1]}]"))
                
            
            print("Creating inclusive histograms for trigger", trigger)
            self.histograms[trigger].extend([
                all_rdf.Histo2D(("Inclusive_EtaVsPt_all", "Inclusive_EtaVsPt;|#eta|;p_{T} (GeV);", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "weight"),
                selected_rdf.Histo2D(("Inclusive_EtaVsPt_selected", "Inclusive_EtaVsPt;|#eta_{jet}|;p_{T} (GeV)", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "weight")
            ]
            )
            
            self.histograms[trigger].extend([eta_binned_rdfs[i].Histo1D(("Inclusive_Pt_eta_" + str(eta_bins_for_pt[i][1]), "Inclusive_pT_eta_" + str(eta_bins_for_pt[i][1]) + ";p_{T} (GeV)", bins_pT_n, bins_pT), "Jet_pt", "weight") for i in range(len(eta_binned_rdfs))])
            
        return 0
    
    def do_PFComposition(self):
        for trigger, rdf in self.trigger_rdfs.items():
            all_rdf = rdf
            selected_rdf = (self.__Flag_cut(rdf)).Redefine("Jet_pt", "Jet_pt[Jet_jetId >= 4]").Redefine("Jet_eta", "Jet_eta[Jet_jetId >= 4]")
            
            print("Creating PFComposition histograms for trigger", trigger)
            # TODO: This kind of behaviour of repeating histogram creation could be optimized
            self.histograms[trigger].extend([
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfilePt_all", "PFComposition_EtaVsPtVsProfilePt;|#eta|;p_{T} (GeV);p_{T} (GeV);", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_pt", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileRho_all", "PFComposition_EtaVsPtVsProfileRho;|#eta|;p_{T} (GeV);#rho;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Rho_fixedGridRhoFastjetAll", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNHF_all", "PFComposition_EtaVsPtVsProfileNHF;|#eta|;p_{T} (GeV);NHF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_neHEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNEF_all", "PFComposition_EtaVsPtVsProfileNEF;|#eta|;p_{T} (GeV);NEF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_neEmEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCHF_all", "PFComposition_EtaVsPtVsProfileCHF;|#eta|;p_{T} (GeV);CHF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_chHEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCEF_all", "PFComposition_EtaVsPtVsProfileCEF;|#eta|;p_{T} (GeV);CEF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_chEmEF", "weight"),
                all_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileMUF_all", "PFComposition_EtaVsPtVsProfileMUF;|#eta|;p_{T} (GeV);MUF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_muEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfilePt_selected", "PFComposition_EtaVsPtVsProfilePt|#eta_{jet}|;p_{T} (GeV);p_{T} (GeV);", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_pt", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileRho_selected", "PFComposition_EtaVsPtVsProfileRho|#eta_{jet}|;p_{T} (GeV);#rho;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Rho_fixedGridRhoFastjetAll", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNHF_selected", "PFComposition_EtaVsPtVsProfileNHF|#eta_{jet}|;p_{T} (GeV);NHF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_neHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileNEF_selected", "PFComposition_EtaVsPtVsProfileNEF|#eta_{jet}|;p_{T} (GeV);NEF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_neEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCHF_selected", "PFComposition_EtaVsPtVsProfileCHF|#eta_{jet}|;p_{T} (GeV);CHF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_chHEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileCEF_selected", "PFComposition_EtaVsPtVsProfileCEF|#eta_{jet}|;p_{T} (GeV);CEF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_chEmEF", "weight"),
                selected_rdf.Profile2D(("PFComposition_EtaVsPtVsProfileMUF_selected", "PFComposition_EtaVsPtVsProfileMUF|#eta_{jet}|;p_{T} (GeV);MUF;", bins_eta_n, bins_eta, bins_pT_n, bins_pT), "Jet_eta", "Jet_pt", "Jet_muEF", "weight")
            ])
            
        return 0
    
    def do_DB(self, system : str = "dijet"):
        for trigger, rdf in self.trigger_rdfs.items():
            if system == "dijet":
                db_rdf = self.__Flag_cut(self.__dijet_cut(rdf))
                pT_binLabels = ["average_Pt_dijet", "Jet_pt_tag", "Jet_pt_probe"]
            elif system == "multijet":
                db_rdf = self.__Flag_cut(self.__multijet_cut(rdf))
                pT_binLabels = ["average_Pt_multijet", "Jet_pt_lead", "pt_recoil"]
            else:
                raise ValueError("System not recognized")

            db_rdf = (db_rdf.Define(f"response_DB_{system}", f"(1.0 + Asymmetry_{system}) / (1.0 - Asymmetry_{system})")
                    )
            
            print("Creating DB histograms for trigger", trigger)
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"DB_{system}_Response", "DB_"+ str(system) + "_response;response;N_{events}", 100, 0.0, 2.0), f"response_DB_{system}", "weight"),
                db_rdf.Histo2D((f"DB_{system}_EtaVsResponse", "DB_"+ str(system) + "_EtaVsResponse;|#eta|;response", bins_eta_n, bins_eta, 100, 0.0, 2.0), f"Jet_eta_tag_{system}", f"response_DB_{system}", "weight")
            ])
            for label in pT_binLabels:
                self.histograms[trigger].append(
                    db_rdf.Histo2D((f"DB_{system}_PtVsResponse_{label}", "DB_"+str(system)+"_PtVsResponse_"+label+";p_{T, "+label+"} (GeV);response;N_{events}",
                                    bins_pT_n, bins_pT, 100, 0.0, 2.0),
                                    label, f"response_DB_{system}", "weight"),
                )

        return 0
    
    def do_MPF(self, system : str = "dijet"):
        for trigger, rdf in self.trigger_rdfs.items():
            if system == "dijet":
                db_rdf = self.__Flag_cut(self.__dijet_cut(rdf))
                pT_binLabels = ["average_Pt_dijet", "Jet_pt_tag", "Jet_pt_probe"]
            elif system == "multijet":
                db_rdf = self.__Flag_cut(self.__multijet_cut(rdf))
                pT_binLabels = ["average_Pt_multijet", "Jet_pt_lead", "pt_recoil"]
            else:
                raise ValueError("System not recognized")

            db_rdf = (db_rdf.Define(f"response_MPF_{system}", f"(1.0 + B_{system}) / (1.0 - B_{system})")
                    )
            
            print("Creating MPF histograms for trigger", trigger)
            self.histograms[trigger].extend([
                db_rdf.Histo1D((f"MPF_{system}_Response", "MPF_"+ str(system) + "_response;response;N_{events}", 100, 0.0, 2.0), f"response_MPF_{system}", "weight"),
                db_rdf.Histo2D((f"DB_{system}_EtaVsResponse", "DB_"+ str(system) + "_EtaVsResponse;|#eta|;response", bins_eta_n, bins_eta, 100, 0.0, 2.0), f"Jet_eta_tag_{system}", f"response_MPF_{system}", "weight")
            ])
            for label in pT_binLabels:
                self.histograms[trigger].append(
                    db_rdf.Histo2D((f"MPF_{system}_PtVsResponse_{label}", "MPF_"+str(system)+"_PtVsResponse_"+label+";p_{T, "+label+"} (GeV);response;N_{events}",
                                    bins_pT_n, bins_pT, 100, 0.0, 2.0),
                                    label, f"response_MPF_{system}", "weight"),
                )

        return 0
    
    
    def __dijet_cut(self, rdf : RDataFrame):
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
    
    def __multijet_cut(self, rdf : RDataFrame):
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
                    .Define("sortedArgs", "ROOT::VecOps::Argsort(Jet_pt)")
                    .Redefine("leading_idx", "sortedArgs[Jet_pt.size() - 1]")
                    .Define("second_idx", "sortedArgs[Jet_pt.size() - 2]")
                    .Define("third_idx", "sortedArgs[Jet_pt.size() - 3]")
                    .Define("Jet_pt_lead", "Jet_pt[leading_idx]")
                    .Define("pt_recoil", "ROOT::VecOps::Sum(Jet_pt) - Jet_pt_lead")
                    .Define("eta_recoil", "ROOT::VecOps::Sum(Jet_eta) - Jet_eta[leading_idx]")
                    .Define("phi_recoil", "ROOT::VecOps::Sum(Jet_phi) - Jet_phi[leading_idx]")
                    .Define("mass_recoil", "ROOT::VecOps::Sum(Jet_mass) - Jet_mass[leading_idx]")
                    .Define("average_Pt_multijet", "(Jet_pt_lead + pt_recoil) / (float(Jet_pt.size()))")
                    .Define("Jet_eta_tag_multijet", "eta_recoil") # Also, scuffed
                    # Is this deltaPhi good? Phi_recoil sussy
                    .Define("deltaPhi_multijet", "abs(Jet_phi[leading_idx] - phi_recoil) < 3.141592653 ? abs(Jet_phi[leading_idx] - phi_recoil) : 2 * 3.141592653 - abs(Jet_phi[leading_idx] - phi_recoil)")
                    # Filter the lead and recoil system
                    .Filter(f"Jet_pt[leading_idx] > {min_pt} && abs(Jet_eta[leading_idx]) < {lead_eta} && Jet_jetId[leading_idx] >= 4")
                    .Filter(f"Jet_pt[second_idx] > {min_pt} && abs(Jet_eta[second_idx]) < {recoil_eta} && Jet_jetId[second_idx] >= 4")
                    .Filter(f"Jet_pt[third_idx] > {min_pt} && abs(Jet_eta[third_idx]) < {recoil_eta} && Jet_jetId[third_idx] >= 4")
                    .Filter(f"Jet_pt[second_idx] < 0.6 * pt_recoil")
                    .Filter(f"Jet_pt[third_idx] < 0.6 * pt_recoil")
                    .Filter(f"abs(deltaPhi_multijet - 3.141592653) < {delta_phi}")
                    # Asymmetry of the system
                    .Define("Asymmetry_multijet", "(Jet_pt_lead - pt_recoil) / (Jet_pt_lead + pt_recoil)")
                    .Define("B_multijet", "RawPuppiMET_pt * cos(ROOT::VecOps::DeltaPhi(Jet_phi[leading_idx], RawPuppiMET_phi)) / (Jet_pt_lead + pt_recoil)")
                    # TODO: Add jet veto
                    )
        self.histograms["all"].extend([rdf_dijet.Histo1D(("Asymmetry", "Asymmetry;Asymmetry;N_{events}", 100, -1.0, 1.0), f"Asymmetry_multijet", "weight"),
                                       rdf_dijet.Histo1D(("pt_recoil", "pt_recoil;p_{T, recoil} (GeV);N_{events}", 100, 0.0, 1000.0), "pt_recoil", "weight"),
                                       rdf_dijet.Histo1D(("Jet_pt_lead", "Jet_pt_lead;p_{T, lead} (GeV);N_{events}", 100, 0.0, 1000.0), "Jet_pt_lead", "weight"),
                                       rdf_dijet.Histo1D(("deltaPhi_multijet", "deltaPhi_multijet;#Delta#phi_{lead, recoil};N_{events}", 100, 0.0, 3.141592653), "deltaPhi_multijet", "weight"),
                                       rdf_dijet.Histo1D(("average_Pt_multijet", "average_Pt_multijet;average p_{T} (GeV);N_{events}", 100, 0.0, 1000.0), "average_Pt_multijet", "weight"),])

        return rdf_dijet
        