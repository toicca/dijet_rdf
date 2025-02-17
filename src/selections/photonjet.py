import ROOT

def init_photonjet(rdf, jet_columns):
    ROOT.gInterpreter.Declare("""
    #ifndef EGAMMA_IDXS
    #define EGAMMA_IDXS
    
    int findPhotonIdx(ROOT::RVec<float> Photon_eta, ROOT::RVec<float> Photon_pt,
                        ROOT::RVec<int> Photon_cutBased, ROOT::RVec<float> Photon_hoe,
                        ROOT::RVec<float> Photon_r9) {
        for (int i = 0; i < Photon_pt.size(); i++) {
            if (abs(Photon_eta[i]) < 1.3 && Photon_pt[i] > 15 && Photon_cutBased[i] == 3 &&
                Photon_hoe[i] < 0.02148 && Photon_r9[i] > 0.94 && Photon_r9[i] < 1.00) {
                return i;
            }
        }
                                
        return -1;
    }
                                
    std::pair<int, int> findJetIdx(ROOT::RVec<float> Jet_eta, ROOT::RVec<float> Jet_pt,
                        ROOT::RVec<float> Jet_phi, ROOT::RVec<int> Jet_jetId,
                        int Photon_jetIdx, float Photon_phi) {
        int idx1 = -1;
        int idx2 = -1;
                                
        // Find the probe jet
        for (int i = 0; i < Jet_pt.size(); i++) {
            if (Jet_pt[i] > 12 && Jet_jetId[i] >= 4
                    && i != Photon_jetIdx
                    && abs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Photon_phi)) > 2.7) {
                idx1 = i;
                break;
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

    rdf = (rdf.Define("Tag_idx_temp", "findPhotonIdx(Photon_eta, Photon_pt, Photon_cutBased, \
            Photon_hoe, Photon_r9)")
            .Filter("Tag_idx_temp >= 0", "Photon found")
            .Define("Jet_indices_temp", "findJetIdx(Jet_eta, Jet_pt, Jet_phi, Jet_jetId, \
                    Photon_jetIdx[Tag_idx_temp], Photon_phi[Tag_idx_temp])")
            .Define("Probe_idx_temp", "Jet_indices_temp.first")
            .Filter("Jet_vetoed[Probe_idx_temp] == 0", "Probe not vetoed")
            .Define("Activity_idx_temp", "Jet_indices_temp.second")
            .Filter("Probe_idx_temp >= 0", "Jet found")
            .Define("Tag_pt", f"Photon_pt[Tag_idx_temp]")
            .Define("Tag_rawPt", "Tag_pt")
            .Define("Tag_eta", f"Photon_eta[Tag_idx_temp]")
            .Define("Tag_phi", f"Photon_phi[Tag_idx_temp]")
            .Define("Tag_mass", "0.0")
            .Define("Tag_label", "2")
            .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    for column in jet_columns:
        rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    return rdf
