import ROOT

def init_zjet(rdf, jet_columns):

    rdf = (rdf.Define("goodMuon_pt", "Muon_pt[Muon_tightId >= 1 && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_pfRelIso03_all < 0.15]")
                .Define("goodMuon_eta", "Muon_eta[Muon_tightId >= 1 && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_pfRelIso03_all < 0.15]")
                .Define("goodMuon_phi", "Muon_phi[Muon_tightId >= 1 && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_pfRelIso03_all < 0.15]")
                .Define("goodMuon_mass", "Muon_mass[Muon_tightId >= 1 && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_pfRelIso03_all < 0.15]")
                .Define("goodMuon_charge", "Muon_charge[Muon_tightId >= 1 && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_pfRelIso03_all < 0.15]")
                .Define("goodElectron_pt", "Electron_pt[Electron_cutBased >= 3 && Electron_pt > 15 && abs(Electron_eta) < 2.4 && Electron_pfRelIso03_all < 0.15]")
                .Define("goodElectron_eta", "Electron_eta[Electron_cutBased >= 3 && Electron_pt > 15 && abs(Electron_eta) < 2.4 && Electron_pfRelIso03_all < 0.15]")
                .Define("goodElectron_phi", "Electron_phi[Electron_cutBased >= 3 && Electron_pt > 15 && abs(Electron_eta) < 2.4 && Electron_pfRelIso03_all < 0.15]")
                .Define("goodElectron_mass", "Electron_mass[Electron_cutBased >= 3 && Electron_pt > 15 && abs(Electron_eta) < 2.4 && Electron_pfRelIso03_all < 0.15]")
                .Define("goodElectron_charge", "Electron_charge[Electron_cutBased >= 3 && Electron_pt > 15 && abs(Electron_eta) < 2.4 && Electron_pfRelIso03_all < 0.15]")
                .Filter("goodMuon_pt.size() >= 2 ? goodElectron_pt.size() == 0 : (goodElectron_pt.size() >= 2 && goodMuon_pt.size() == 0)", "Only one lepton flavor")
    )
        
    ROOT.gInterpreter.Declare("""
    #ifndef ZJET_IDXS
    #define ZJET_IDXS
                                
    std::pair<int, int> findMuonIdxs(ROOT::RVec<float> Muon_eta, ROOT::RVec<float> Muon_pt,
                        ROOT::RVec<float> Muon_pfRelIso03_all, ROOT::RVec<int> Muon_tightId,
                        ROOT::RVec<float> Muon_charge) {
        int idx1 = -1;
        int idx2 = -1;
        for (int i = 0; i < Muon_pt.size(); i++) {
            if (abs(Muon_eta[i]) < 2.4 && Muon_pfRelIso03_all[i] < 0.15 &&
                Muon_tightId[i]) {
                // Leading muon pt>20, subleading pt>10
                if (idx1 == -1 && Muon_pt[i] > 20) {
                    idx1 = i;
                } else if (idx2 == -1 && Muon_charge[i] != Muon_charge[idx1] &&
                        Muon_pt[i] > 10) {
                    idx2 = i;
                    break;
                }
            }
        }
                                
        return std::make_pair(idx1, idx2);
    }

    std::pair<int, int> findJetIdx(ROOT::RVec<float> Jet_eta, ROOT::RVec<float> Jet_pt,
                        ROOT::RVec<float> Jet_phi, ROOT::RVec<int> Jet_jetId,
                        float Z_eta, float Z_phi) {

        if (Jet_pt[0] < 12 || Jet_jetId[0] < 4 || abs(ROOT::VecOps::DeltaPhi(Jet_phi[0], Z_phi)) < 2.7) {
            return std::make_pair(-1, -1);
        }
                              
        if (Jet_pt[1] < 12 || Jet_jetId[1] < 4) {
            return std::make_pair(0, -1);
        }
                              
        return std::make_pair(0, 1);
    }
                                
    #endif
    """)
    rdf = (rdf.Filter("nMuon > 1", "nMuon > 1")
            .Define("Muon_idx_temp", "findMuonIdxs(Muon_eta, Muon_pt, Muon_pfRelIso03_all, \
                    Muon_tightId, Muon_charge)")
            .Filter("Muon_idx_temp.first >= 0 && Muon_idx_temp.second >= 0", "Two muons found")
            .Define("Z_4vec_temp",
                "ROOT::Math::PtEtaPhiMVector(Muon_pt[Muon_idx_temp.first], \
                        Muon_eta[Muon_idx_temp.first], Muon_phi[Muon_idx_temp.first], \
                        Muon_mass[Muon_idx_temp.first]) + \
                        ROOT::Math::PtEtaPhiMVector(Muon_pt[Muon_idx_temp.second], \
                        Muon_eta[Muon_idx_temp.second], Muon_phi[Muon_idx_temp.second], \
                        Muon_mass[Muon_idx_temp.second])")
            .Define("Tag_pt", "static_cast<float>(Z_4vec_temp.Pt())")
            .Define("Tag_rawPt", "Tag_pt")
            .Define("Tag_eta", "static_cast<float>(Z_4vec_temp.Eta())")
            .Define("Tag_phi", "static_cast<float>(Z_4vec_temp.Phi())")
            .Define("Tag_mass", "static_cast<float>(Z_4vec_temp.M())")
            .Define("Tag_label", "1")
            .Define("Jet_indices_temp", "findJetIdx(Jet_eta, Jet_pt, Jet_phi, Jet_jetId, \
                    Tag_eta, Tag_phi)")
            .Define("Probe_idx_temp", "Jet_indices_temp.first")
            .Filter("Jet_vetoed[Probe_idx_temp] == 0", "Jet not vetoed")
            .Define("Activity_idx_temp", "Jet_indices_temp.second")
            .Filter("Probe_idx_temp >= 0", "Jet found")
            .Filter("Tag_pt > 12", "Z pT > 12")
            .Filter("Tag_mass > 71.1876 && Tag_mass < 111.1876", "Z mass window")
            .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    for column in jet_columns:
        rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    return rdf