#ifndef MULTIJET_IDXS
#define MULTIJET_IDXS

#include "ROOT/RVec.hxx"

ROOT::RVec<int> findRecoilJetIdxs(const ROOT::RVec<float>& Jet_pt,
                                  const ROOT::RVec<float>& Jet_eta,
                                  const ROOT::RVec<float>& Jet_phi,
                                  const ROOT::RVec<float>& Jet_mass,
                                  const ROOT::RVec<int>& Jet_jetId);

#endif
