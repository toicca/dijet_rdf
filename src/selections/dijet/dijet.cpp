#include "dijet.h"
#include <random>
#include <cmath>
#include "ROOT/VecOps.hxx"

std::pair<std::pair<int, int>, int> findTagProbeIdxs(
    ROOT::RVec<float> Jet_eta,
    ROOT::RVec<float> Jet_pt,
    ROOT::RVec<float> Jet_phi,
    ROOT::RVec<int> Jet_jetId
) {
    // Random number generator setup
    static std::random_device rd;  // Non-deterministic random source
    static std::mt19937 gen(rd()); // Seeded Mersenne Twister PRNG
    static std::uniform_int_distribution<int> dist(0, 1); // Uniform distribution: 0 or 1

    int idx1 = dist(gen); // Generate random 0 or 1
    int idx2 = 1 - idx1;

    // Check that the tag is in barrel
    if (abs(Jet_eta[idx1]) > 1.3 || Jet_pt[idx1] < 12 || Jet_jetId[idx1] < 4) {
        return std::make_pair(std::make_pair(-1, -1), -1);
    }

    // Check that the probe is back-to-back with the tag
    if (abs(ROOT::VecOps::DeltaPhi(Jet_phi[idx2], Jet_phi[idx1])) < 2.7 ||
        Jet_pt[idx2] < 12 || Jet_jetId[idx2] < 4) {
        return std::make_pair(std::make_pair(-1, -1), -1);
    }

    // Find the activity jet
    int idx3 = -1;
    for (int i = 0; i < Jet_pt.size(); i++) {
        if (i != idx1 && i != idx2 && Jet_pt[i] > 12 && Jet_jetId[i] >= 4) {
            idx3 = i;
            break;
        }
    }

    return std::make_pair(std::make_pair(idx1, idx2), idx3);
}