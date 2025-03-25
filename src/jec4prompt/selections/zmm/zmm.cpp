#include "zmm.h"
#include "ROOT/RVec.hxx"
#include <cmath>
#include "Math/Vector4D.h"

std::pair<int, int> findMuonIdxs(const ROOT::RVec<float>& Muon_pt,
                                 const ROOT::RVec<float>& Muon_eta,
                                 const ROOT::RVec<float>& Muon_phi,
                                 const ROOT::RVec<float>& Muon_mass,
                                 const ROOT::RVec<float>& Muon_charge) {
    int idx1 = -1;
    int idx2 = -1;

    float Zmass = 91.1876;
    float mtemp = 0.;

    for (int i = 0; i < Muon_pt.size(); i++) {
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

ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Muon_eta,
                           const ROOT::RVec<float>& Muon_phi,
                           const ROOT::RVec<float>& trg_eta,
                           const ROOT::RVec<float>& trg_phi,
                           const ROOT::RVec<int>& trg_filterBits,
                           const ROOT::RVec<int>& trg_id) {
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

std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_pt,
                                const ROOT::RVec<float>& Jet_eta,
                                const ROOT::RVec<float>& Jet_phi,
                                const ROOT::RVec<float>& Muon_eta,
                                const ROOT::RVec<float>& Muon_phi) {
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
            } else if (idx2 == -1 && Jet_pt[i] > 15) {
                idx2 = i;
                break;
            }
        }
    }
    return std::make_pair(idx1, idx2);
}

ROOT::RVec<bool> passMuMET(const ROOT::RVec<float>& Muon_eta,
                        const ROOT::RVec<float>& Muon_phi,
                        const ROOT::RVec<float>& Jet_pt,
                        const ROOT::RVec<float>& Jet_eta,
                        const ROOT::RVec<float>& Jet_phi,
                        const ROOT::RVec<float>& Jet_rawFactor,
                        const ROOT::RVec<float>& Jet_muonSubtrFactor) {
    ROOT::RVec<bool> passMuMET(Jet_eta.size(), true);

    // dR check
    for (int i = 0; i < Jet_eta.size(); i++) {
        for (int j = 0; j < Muon_eta.size(); j++) {
            if (ROOT::VecOps::DeltaR(Jet_eta[i], Muon_eta[j], Jet_phi[i], Muon_phi[j]) < 0.3) {
                passMuMET[i] = false;
                break;
            }
        }
    }

    // pT check
    // Require that the corrected pT of the jet is greater than 15 GeV after muon subtraction
    // TODO: Should the muon be subtracted from the probe jet?
    for (int i = 0; i < Jet_eta.size(); i++) {
        if (!passMuMET[i]) continue;

        float rpt = (1.0 - Jet_rawFactor[i]) * Jet_pt[i];
        float musubtrpt = rpt * (1.0 - Jet_muonSubtrFactor[i]);
        float mupt = rpt - musubtrpt;

        if (Jet_pt[i] - mupt < 15) {
            passMuMET[i] = false;
        }
    }


    return passMuMET;
}