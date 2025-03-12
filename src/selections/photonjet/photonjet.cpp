#include "photonjet.h"
#include "ROOT/RVec.hxx"
#include <cmath>

using namespace ROOT;

ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Photon_eta,
                           const ROOT::RVec<float>& Photon_phi,
                           const ROOT::RVec<float>& trg_eta,
                           const ROOT::RVec<float>& trg_phi,
                           const ROOT::RVec<int>& trg_id) {
    ROOT::RVec<bool> trgObj(Photon_eta.size(), false);
    for (int i = 0; i < trgObj.size(); i++) {
        for (int j = 0; j < trg_eta.size(); j++) {
            if (ROOT::VecOps::DeltaR(Photon_eta[i], trg_eta[j], Photon_phi[i], trg_phi[j]) < 0.3) {
                if (trg_id[j] == 22) {
                    trgObj[i] = true;
                    break;
                }
            }
        }
    }
    return trgObj;
}

std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_eta,
                                const ROOT::RVec<float>& Jet_phi,
                                const ROOT::RVec<float>& Photon_eta,
                                const ROOT::RVec<float>& Photon_phi) {
    int idx1 = -1;
    int idx2 = -1;
    for (int i = 0; i < Jet_eta.size(); i++) {
        bool badJet = false;
        for (int j = 0; j < Photon_eta.size(); j++) {
            if (ROOT::VecOps::DeltaR(Jet_eta[i], Photon_eta[j], Jet_phi[i], Photon_phi[j]) < 0.3) {
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