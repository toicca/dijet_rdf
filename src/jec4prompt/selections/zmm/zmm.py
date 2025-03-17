import ROOT

def init_zmm(rdf, jet_columns, state):
    path = state.module_dir
    path = path / "selections" / "zmm"
    cpp_path = path / "zmm.cpp"
    so_path = path / "zmm_cpp.so"
    h_path = path / "zmm.h"

    # Compile and load the C++ code
    ROOT.gInterpreter.ProcessLine(f'.L {cpp_path}+')
    ROOT.gSystem.Load(str(so_path))
    ROOT.gInterpreter.Declare(f'#include "{h_path}"')

    # Good muon selection
    rdf = (rdf.Define("muonMask", "Muon_pt > 8 && Muon_pfIsoId >= 4 && Muon_pfRelIso04_all < 0.15 && Muon_tightId")
            .Define("goodMuon_pt", "Muon_pt[muonMask]")
            .Define("goodMuon_eta", "Muon_eta[muonMask]")
            .Define("goodMuon_charge", "Muon_charge[muonMask]")
            .Define("goodMuon_phi", "Muon_phi[muonMask]")
            .Define("goodMuon_mass", "Muon_mass[muonMask]")
    )

    # Trigger selected muons
    rdf = (rdf.Define("trigMask", "hasTrgObj(goodMuon_eta, goodMuon_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)")
            .Define("selMuon_pt", "goodMuon_pt[trigMask]")
            .Define("selMuon_eta", "goodMuon_eta[trigMask]")
            .Define("selMuon_charge", "goodMuon_charge[trigMask]")
            .Define("selMuon_phi", "goodMuon_phi[trigMask]")
            .Define("selMuon_mass", "goodMuon_mass[trigMask]")
            .Filter("selMuon_pt.size() > 1 && selMuon_pt.size() < 4",
                    "2-3 tight muons with trigger match")
            .Define("Muon_idx_temp", "findMuonIdxs(selMuon_pt, selMuon_eta, selMuon_phi, selMuon_mass, selMuon_charge)")
            .Filter("Muon_idx_temp.first >= 0 && Muon_idx_temp.second >= 0", 
                    "Two muons found")
            .Filter("selMuon_pt[Muon_idx_temp.first] > 20 && selMuon_pt[Muon_idx_temp.second] > 10",
                    "Leading muon pT > 20, subleading muon pT > 10")
            .Filter("fabs(selMuon_eta[Muon_idx_temp.first]) <= 2.3 && fabs(selMuon_eta[Muon_idx_temp.second]) <= 2.3",
                    "|Muon eta| <= 2.3")
    )

    # Z boson reconstruction / Tag object
    rdf = (rdf.Define("Z_4vec_temp",
                "ROOT::Math::PtEtaPhiMVector(selMuon_pt[Muon_idx_temp.first], \
                        selMuon_eta[Muon_idx_temp.first], selMuon_phi[Muon_idx_temp.first], \
                        selMuon_mass[Muon_idx_temp.first]) + \
                        ROOT::Math::PtEtaPhiMVector(selMuon_pt[Muon_idx_temp.second], \
                        selMuon_eta[Muon_idx_temp.second], selMuon_phi[Muon_idx_temp.second], \
                        selMuon_mass[Muon_idx_temp.second])")
            .Define("Tag_pt", "static_cast<float>(Z_4vec_temp.Pt())")
            .Define("Tag_rawPt", "Tag_pt")
            .Define("Tag_eta", "static_cast<float>(Z_4vec_temp.Eta())")
            .Define("Tag_phi", "static_cast<float>(Z_4vec_temp.Phi())")
            .Define("Tag_mass", "static_cast<float>(Z_4vec_temp.M())")
            .Define("Tag_label", "1")
            .Filter("Tag_pt > 12",
                    "Z pT > 12")
            .Filter("Tag_mass > 71.1876 && Tag_mass < 111.1876",
                    "Z mass window, +-20 GeV")
    )

    # Probe jet selection
    rdf = (rdf.Define("JetMuon_idx_temp", "findJetIdxs(Jet_pt, Jet_eta, Jet_phi, goodMuon_eta, goodMuon_phi)")
            .Define("Probe_idx_temp", "JetMuon_idx_temp.first")
            .Filter("Probe_idx_temp >= 0",
                    "Jet found")
            .Filter("Jet_pt[Probe_idx_temp] > 12",
                    "Leading jet pT > 12")
            # .Filter("(1.0-Jet_rawFactor[Probe_idx_temp])*Jet_pt[Probe_idx_temp] > 12",
                    # "Leading jet raw pT > 12")
            .Filter("Jet_jetId[Probe_idx_temp] >= 4",
                    "Leading jet Id >= 4")
            .Filter("Jet_vetoed[Probe_idx_temp] == 0",
                    "Jet not vetoed")
            .Filter("fabs(ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)) > 2.7",
                    "|dPhi(Z,jet)| > 2.7")
            .Define("Activity_idx_temp", "JetMuon_idx_temp.second")
            .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    rdf = rdf.Define("Activity_denom", "Tag_pt")

    for column in jet_columns:
        rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    return rdf