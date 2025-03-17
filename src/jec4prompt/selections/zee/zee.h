#ifndef ZJET_IDXS
#define ZJET_IDXS

#include "ROOT/RVec.hxx"
#include <utility>

std::pair<int, int> findElectronIdxs(const ROOT::RVec<float>& Electron_pt,
                                     const ROOT::RVec<float>& Electron_eta,
                                     const ROOT::RVec<float>& Electron_phi,
                                     const ROOT::RVec<float>& Electron_mass,
                                     const ROOT::RVec<float>& Electron_charge);

std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_eta,
                                const ROOT::RVec<float>& Jet_phi,
                                const ROOT::RVec<float>& Electron_eta,
                                const ROOT::RVec<float>& Electron_phi);

ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Electron_eta,
                           const ROOT::RVec<float>& Electron_phi,
                           const ROOT::RVec<float>& trg_eta,
                           const ROOT::RVec<float>& trg_phi,
                           const ROOT::RVec<int>& trg_filterBits,
                           const ROOT::RVec<int>& trg_id);

#endif
