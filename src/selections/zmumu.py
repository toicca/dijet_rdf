import ROOT

def init_zmumu(rdf, jet_columns):
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
        int idx1 = -1;
        int idx2 = -1;
        
        // Find the probe jet
        for (int i = 0; i < Jet_pt.size(); i++) {
            if (Jet_pt[i] > 12 && Jet_jetId[i] >= 4 && abs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Z_phi)) > 2.7) {
                if (idx1 == -1) {
                    idx1 = i;
                    break;
                }
            }
        }

        // Find the activity jet
        for (int i = 0; i < Jet_pt.size(); i++) {
            if (i != idx1 && Jet_pt[i] > 12 && Jet_jetId[i] >= 4) {
                idx2 = i;
                break;
            }
        }
                                
        return std::make_pair(idx1, idx2);
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
            .Define("Activity_idx_temp", "Jet_indices_temp.second")
            .Filter("Probe_idx_temp >= 0", "Jet found")
            .Filter("Tag_pt > 12", "Z pT > 12")
            .Filter("Tag_mass > 71.1876 && Tag_mass < 111.1876", "Z mass window")
    )

    for column in jet_columns:
        rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    return rdf