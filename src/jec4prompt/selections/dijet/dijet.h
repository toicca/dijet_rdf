#ifndef DIJET_H
#define DIJET_H

#include <utility>
#include "ROOT/RVec.hxx"

std::pair<std::pair<int, int>, int> findTagProbeIdxs(
    ROOT::RVec<float> Jet_eta,
    ROOT::RVec<float> Jet_pt,
    ROOT::RVec<float> Jet_phi,
    ROOT::RVec<int> Jet_jetId
);

#endif
