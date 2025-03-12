#include "multijet.h"
#include "ROOT/VecOps.hxx"
#include <cmath>

ROOT::RVec<int> findRecoilJetIdxs(const ROOT::RVec<float>& Jet_pt,
                                  const ROOT::RVec<float>& Jet_eta,
                                  const ROOT::RVec<float>& Jet_phi,
                                  const ROOT::RVec<float>& Jet_mass,
                                  const ROOT::RVec<int>& Jet_jetId) {
    ROOT::RVec<int> idxs;

    for (int i = 1; i < Jet_pt.size(); i++) {
        if (Jet_pt[i] > 30 && abs(Jet_eta[i]) < 2.5 && Jet_jetId[i] >= 4
            && abs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Jet_phi[0])) > 1.0) {
            idxs.push_back(i);
        }
    }

    return idxs;
}