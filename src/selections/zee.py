import ROOT

def init_zee(rdf, jet_columns):
    ROOT.gInterpreter.Declare("""
    #ifndef ZJET_IDXS
    #define ZJET_IDXS
                                
    std::pair<int, int> findElectronIdxs(ROOT::RVec<float> Electron_eta, ROOT::RVec<float> Electron_pt,
                        ROOT::RVec<float> Electron_pfRelIso03_all, ROOT::RVec<int> Electron_cutBased,
                        ROOT::RVec<float> Electron_charge) {
        int idx1 = -1;
        int idx2 = -1;
        for (int i = 0; i < Electron_pt.size(); i++) {
            if (abs(Electron_eta[i]) < 2.4 && Electron_pfRelIso03_all[i] < 0.15 &&
                Electron_cutBased[i] >= 4) {
                if (idx1 == -1 && Electron_pt[i] > 15) {
                    idx1 = i;
                } else if (idx2 == -1 && Electron_charge[i] != Electron_charge[idx1] &&
                        Electron_pt[i] > 15) {
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
    rdf = (rdf.Filter("nElectron > 1", "nElectron > 1")
            .Define("Electron_idx_temp", "findElectronIdxs(Electron_eta, Electron_pt, Electron_pfRelIso03_all, \
                    Electron_cutBased, Electron_charge)")
            .Filter("Electron_idx_temp.first >= 0 && Electron_idx_temp.second >= 0", "Two electrons found")
            .Define("Z_4vec_temp",
                "ROOT::Math::PtEtaPhiMVector(Electron_pt[Electron_idx_temp.first], \
                        Electron_eta[Electron_idx_temp.first], Electron_phi[Electron_idx_temp.first], \
                        Electron_mass[Electron_idx_temp.first]) + \
                        ROOT::Math::PtEtaPhiMVector(Electron_pt[Electron_idx_temp.second], \
                        Electron_eta[Electron_idx_temp.second], Electron_phi[Electron_idx_temp.second], \
                        Electron_mass[Electron_idx_temp.second])")
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