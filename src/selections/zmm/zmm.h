#ifndef ZJET_IDXS
#define ZJET_IDXS

#include "ROOT/RVec.hxx"
#include <utility>

std::pair<int, int> findMuonIdxs(const ROOT::RVec<float>& Muon_pt,
                                 const ROOT::RVec<float>& Muon_eta,
                                 const ROOT::RVec<float>& Muon_phi,
                                 const ROOT::RVec<float>& Muon_mass,
                                 const ROOT::RVec<float>& Muon_charge);

ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Muon_eta,
                           const ROOT::RVec<float>& Muon_phi,
                           const ROOT::RVec<float>& trg_eta,
                           const ROOT::RVec<float>& trg_phi,
                           const ROOT::RVec<int>& trg_filterBits,
                           const ROOT::RVec<int>& trg_id);

std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_eta,
                                const ROOT::RVec<float>& Jet_phi,
                                const ROOT::RVec<float>& Muon_eta,
                                const ROOT::RVec<float>& Muon_phi);

#endif
