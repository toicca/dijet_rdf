#include "multijet.h"
#include "ROOT/RVec.hxx"
#include <cmath>

ROOT::RVec<int> findRecoilJetIdxs(const ROOT::RVec<float>& Jet_pt,
                                  const ROOT::RVec<float>& Jet_eta,
                                  const ROOT::RVec<float>& Jet_phi,
                                  const ROOT::RVec<float>& Jet_mass,
                                  const ROOT::RVec<int>& Jet_jetId) {
    ROOT::RVec<int> idxs;

    for (int i = 1; i < Jet_pt.size(); i++) {
        if (Jet_pt[i] > 30 && fabs(Jet_eta[i]) < 2.5
            && fabs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Jet_phi[0])) > 1.0) {
            idxs.push_back(i);
        }
    }

    return idxs;
}

bool multijetVetoForward(const ROOT::RVec<float>& Jet_pt,
                         const ROOT::RVec<float>& Jet_eta) {
    for (int i = 1; i < Jet_eta.size(); i++) {
        if (Jet_pt[i] > 30 && fabs(Jet_eta[i]) >= 2.5) {
            return false;
        }
    }

    return true;
}

bool multijetVetoNear(const ROOT::RVec<float>& Jet_pt,
                      const ROOT::RVec<float>& Jet_eta,
                      const ROOT::RVec<float>& Jet_phi) {
    for (int i = 1; i < Jet_eta.size(); i++) {
        if (Jet_pt[i] > 30 && fabs(Jet_eta[i]) < 2.5
            && fabs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Jet_phi[0])) <= 1.0) {
            return false;
        }
    }

    return true;
}