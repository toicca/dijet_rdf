#include "zee.h"
#include "ROOT/VecOps.hxx"
#include <cmath>
#include "Math/Vector4D.h"

std::pair<int, int> findElectronIdxs(const ROOT::RVec<float>& Electron_pt,
                                     const ROOT::RVec<float>& Electron_eta,
                                     const ROOT::RVec<float>& Electron_phi,
                                     const ROOT::RVec<float>& Electron_mass,
                                     const ROOT::RVec<float>& Electron_charge) {
    int idx1 = -1;
    int idx2 = -1;

    float Zmass = 91.1876;
    float mtemp = 0.;

    for (int i = 0; i < Electron_pt.size(); i++) {
        ROOT::Math::PtEtaPhiMVector ele1(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);
        for (int j = i+1; j < Electron_pt.size(); j++) {
            if (Electron_charge[i] == Electron_charge[j]) continue;
            ROOT::Math::PtEtaPhiMVector ele2(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
            ROOT::Math::PtEtaPhiMVector Z = ele1 + ele2;

            if (abs(mtemp - Zmass) > abs(Z.M() - Zmass)) {
                mtemp = Z.M();
                idx1 = i;
                idx2 = j;
            }
        }
    }

    if (idx1 == -1 || idx2 == -1) return std::make_pair(-1, -1);

    // Sort electrons by pt
    if (Electron_pt[idx1] < Electron_pt[idx2]) {
        int temp = idx1;
        idx1 = idx2;
        idx2 = temp;
    }

    return std::make_pair(idx1, idx2);
}

std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_eta,
                                const ROOT::RVec<float>& Jet_phi,
                                const ROOT::RVec<float>& Electron_eta,
                                const ROOT::RVec<float>& Electron_phi) {
    int idx1 = -1;
    int idx2 = -1;
    for (int i = 0; i < Jet_eta.size(); i++) {
        bool badJet = false;
        for (int j = 0; j < Electron_eta.size(); j++) {
            if (ROOT::VecOps::DeltaR(Jet_eta[i], Electron_eta[j], Jet_phi[i], Electron_phi[j]) < 0.3) {
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

ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Electron_eta,
                           const ROOT::RVec<float>& Electron_phi,
                           const ROOT::RVec<float>& trg_eta,
                           const ROOT::RVec<float>& trg_phi,
                           const ROOT::RVec<int>& trg_filterBits,
                           const ROOT::RVec<int>& trg_id) {
    ROOT::RVec<bool> trgObj(Electron_eta.size(), false);
    for (int i = 0; i < Electron_eta.size(); i++) {
        for (int j = 0; j < trg_eta.size(); j++) {
            if (ROOT::VecOps::DeltaR(Electron_eta[i], trg_eta[j], Electron_phi[i], trg_phi[j]) < 0.3) {
                if (trg_id[j] == 11 && (trg_filterBits[j] & 1)) {
                    trgObj[i] = true;
                    break;
                }
            }
        }
    }
    return trgObj;
}