import ROOT

def init_zmm(rdf, jet_columns):
    ROOT.gInterpreter.Declare("""
    #ifndef ZJET_IDXS
    #define ZJET_IDXS

    std::pair<int, int> findMuonIdxs(ROOT::RVec<float> Muon_pt, ROOT::RVec<float> Muon_eta,
                        ROOT::RVec<float> Muon_phi, ROOT::RVec<float> Muon_mass,
                        ROOT::RVec<float> Muon_charge) {
        int idx1 = -1;
        int idx2 = -1;

        float Zmass = 91.1876;
        float mtemp = 0.;

        for (int i = 0; i < Muon_pt.size(); i++) {
            // if (Muon_pt[i] < 20) continue;
            // Not having the pt cuts here allows to "fail" based on the z mass?
            ROOT::Math::PtEtaPhiMVector mu1(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
            for (int j = i+1; j < Muon_pt.size(); j++) {
                if (Muon_charge[i] == Muon_charge[j]) continue;
                ROOT::Math::PtEtaPhiMVector mu2(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
                ROOT::Math::PtEtaPhiMVector Z = mu1 + mu2;
                              
                if (abs(mtemp - Zmass) > abs(Z.M() - Zmass)) {
                    mtemp = Z.M();
                    idx1 = i;
                    idx2 = j;
                }
            }
        }

        if (idx1 == -1 || idx2 == -1) return std::make_pair(-1, -1);

        // Sort muons by pt
        if (Muon_pt[idx1] < Muon_pt[idx2]) {
            int temp = idx1;
            idx1 = idx2;
            idx2 = temp;
        }
                                
        return std::make_pair(idx1, idx2);
    }

    ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Muon_eta, const ROOT::RVec<float>& Muon_phi,
        const ROOT::RVec<float>& trg_eta, const ROOT::RVec<float>& trg_phi, const ROOT::RVec<int>& trg_filterBits, const ROOT::RVec<int>& trg_id) {

        ROOT::RVec<bool> trgObj(Muon_eta.size(), false);
        for (int i = 0; i < Muon_eta.size(); i++) {
            for (int j = 0; j < trg_eta.size(); j++) {
                if (ROOT::VecOps::DeltaR(Muon_eta[i], trg_eta[j], Muon_phi[i], trg_phi[j]) < 0.3) {
                    if (trg_id[j] == 13 && (trg_filterBits[j] & (1 << 3))) {
                        trgObj[i] = true;
                        break;
                    }
                }
            }
        }
        return trgObj;
    }


    std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_eta, const ROOT::RVec<float>& Jet_phi,
                              const ROOT::RVec<float>& Muon_eta, const ROOT::RVec<float>& Muon_phi) {
        int idx1 = -1;
        int idx2 = -1;
        for (int i = 0; i < Jet_eta.size(); i++) {
            bool badJet = false;
            for (int j = 0; j < Muon_eta.size(); j++) {
                if (ROOT::VecOps::DeltaR(Jet_eta[i], Muon_eta[j], Jet_phi[i], Muon_phi[j]) < 0.3) {
                    badJet = true;
                    break;
                }
            }
            if (!badJet) {
                if (idx1 == -1) {
                    idx1 = i;
                } else {
                    idx2 = i;
                    break;
                }
            }
        }
        return std::make_pair(idx1, idx2);
    }
                                
    #endif
    """)
    rdf = (rdf.Define("muonMask", "Muon_pt > 8 && Muon_pfIsoId >= 4 && Muon_pfRelIso04_all < 0.15 && Muon_tightId")
            .Define("selMuon_pt", "Muon_pt[muonMask]")
            .Define("selMuon_eta", "Muon_eta[muonMask]")
            .Define("selMuon_charge", "Muon_charge[muonMask]")
            .Define("selMuon_phi", "Muon_phi[muonMask]")
            .Define("selMuon_mass", "Muon_mass[muonMask]")
            .Define("trigMask", "hasTrgObj(selMuon_eta, selMuon_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)")
            .Define("goodMuon_pt", "selMuon_pt[trigMask]")
            .Define("goodMuon_eta", "selMuon_eta[trigMask]")
            .Define("goodMuon_charge", "selMuon_charge[trigMask]")
            .Define("goodMuon_phi", "selMuon_phi[trigMask]")
            .Define("goodMuon_mass", "selMuon_mass[trigMask]")
            .Filter("goodMuon_pt.size() > 1 && goodMuon_pt.size() < 4", "2-3 tight muons with trigger match")
            .Define("Muon_idx_temp", "findMuonIdxs(goodMuon_pt, goodMuon_eta, goodMuon_phi, goodMuon_mass, goodMuon_charge)")
            .Filter("Muon_idx_temp.first >= 0 && Muon_idx_temp.second >= 0", "Two muons found")
            .Filter("goodMuon_pt[Muon_idx_temp.first] > 20 && goodMuon_pt[Muon_idx_temp.second] > 10", "Leading muon pT > 20, subleading muon pT > 12")
            .Filter("abs(goodMuon_eta[Muon_idx_temp.first]) <= 2.3 && abs(goodMuon_eta[Muon_idx_temp.second]) <= 2.3", "Muon eta <= 2.3")
            .Define("Z_4vec_temp",
                "ROOT::Math::PtEtaPhiMVector(goodMuon_pt[Muon_idx_temp.first], \
                        goodMuon_eta[Muon_idx_temp.first], goodMuon_phi[Muon_idx_temp.first], \
                        goodMuon_mass[Muon_idx_temp.first]) + \
                        ROOT::Math::PtEtaPhiMVector(goodMuon_pt[Muon_idx_temp.second], \
                        goodMuon_eta[Muon_idx_temp.second], goodMuon_phi[Muon_idx_temp.second], \
                        goodMuon_mass[Muon_idx_temp.second])")
            .Define("Tag_pt", "static_cast<float>(Z_4vec_temp.Pt())")
            .Define("Tag_rawPt", "Tag_pt")
            .Define("Tag_eta", "static_cast<float>(Z_4vec_temp.Eta())")
            .Define("Tag_phi", "static_cast<float>(Z_4vec_temp.Phi())")
            .Define("Tag_mass", "static_cast<float>(Z_4vec_temp.M())")
            .Define("Tag_label", "1")
            .Define("JetMuon_idx_temp", "findJetIdxs(Jet_eta, Jet_phi, selMuon_eta, selMuon_phi)")
            .Define("Probe_idx_temp", "JetMuon_idx_temp.first")
            .Filter("Tag_pt > 12", "Z pT > 12")
            .Filter("Tag_mass > 71.1876 && Tag_mass < 111.1876", "Z mass window")
            .Filter("Probe_idx_temp >= 0", "Jet found")
            .Filter("Jet_pt[Probe_idx_temp] > 12", "Leading jet pT > 12")
            .Filter("Jet_jetId[Probe_idx_temp] >= 4", "Leading jet Id >= 4")
            .Filter("Jet_vetoed[Probe_idx_temp] == 0", "Jet not vetoed")
            .Filter("abs(ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)) > 2.7", "dPhi(Z,jet) > 2.7")
            .Define("Tag_deltaPhi", "ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)")
            .Define("Activity_idx_temp", "JetMuon_idx_temp.second")
            .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    for column in jet_columns:
        rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    return rdf