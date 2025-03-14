#ifndef EGAMMA_IDXS
#define EGAMMA_IDXS

#include "ROOT/RVec.hxx"
#include <utility>

ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Photon_eta,
                           const ROOT::RVec<float>& Photon_phi,
                           const ROOT::RVec<float>& trg_eta,
                           const ROOT::RVec<float>& trg_phi,
                           const ROOT::RVec<int>& trg_id);

std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_pt,
                                const ROOT::RVec<float>& Jet_eta,
                                const ROOT::RVec<float>& Jet_phi,
                                const ROOT::RVec<float>& Photon_eta,
                                const ROOT::RVec<float>& Photon_phi,
                                const ROOT::RVec<int>& Photon_jetIdx);

ROOT::RVec<int> metJetIdx(const ROOT::RVec<float>& Jet_pt,
                          const ROOT::RVec<float>& Jet_eta,
                          const ROOT::RVec<float>& Jet_phi,
                          float Tag_eta, float Tag_phi);

#endif
