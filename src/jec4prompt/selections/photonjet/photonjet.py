import ROOT


def init_photonjet(rdf, jet_columns, state):
    path = state.module_dir
    path = path / "selections" / "photonjet"
    cpp_path = path / "photonjet.cpp"
    so_path = path / "photonjet_cpp.so"
    h_path = path / "photonjet.h"

    # Compile and load the C++ code
    ROOT.gInterpreter.ProcessLine(f".L {cpp_path}+")
    ROOT.gSystem.Load(str(so_path))
    ROOT.gInterpreter.Declare(f'#include "{h_path}"')

    # Good photons
    rdf = (
        rdf.Define(
            "isGoodPhoton",
            "Photon_cutBased == 3 && \
                      Photon_hoe < 0.02148 && Photon_r9 > 0.94 && Photon_r9 < 1.00",
        )
        .Define("goodPhoton_pt", "Photon_pt[isGoodPhoton]")
        .Define("goodPhoton_eta", "Photon_eta[isGoodPhoton]")
        .Define("goodPhoton_phi", "Photon_phi[isGoodPhoton]")
        .Define("goodPhoton_jetIdx", "Photon_jetIdx[isGoodPhoton]")
        .Define("goodPhoton_seedGain", "Photon_seedGain[isGoodPhoton]")
        .Filter("goodPhoton_pt.size() > 0", "Good photon found")
    )

    # Trigger selected photons
    rdf = (
        rdf.Define(
            "hasTrg_temp",
            "hasTrgObj(goodPhoton_eta, goodPhoton_phi, TrigObj_eta, TrigObj_phi, TrigObj_id)",
        )
        .Define("selPhoton_pt", "goodPhoton_pt[hasTrg_temp]")
        .Define("selPhoton_eta", "goodPhoton_eta[hasTrg_temp]")
        .Define("selPhoton_phi", "goodPhoton_phi[hasTrg_temp]")
        .Define("selPhoton_jetIdx", "goodPhoton_jetIdx[hasTrg_temp]")
        .Define("selPhoton_seedGain", "goodPhoton_seedGain[hasTrg_temp]")
        .Filter("selPhoton_pt.size() > 0", "Trigger selected photon found")
        .Filter("selPhoton_pt[0] > 15", "Trigger selected photon pT > 15")
        .Filter("abs(selPhoton_eta[0]) < 1.3", "Trigger selected photon |eta| < 1.3")
        .Define("Tag_idx_temp", "0")
        .Define("Tag_pt", f"selPhoton_pt[Tag_idx_temp]")
        .Define("Tag_rawPt", "Tag_pt")
        .Redefine(
            "Tag_pt",
            "selPhoton_seedGain[Tag_idx_temp] == 1 ? Tag_pt * 1./1.017: Tag_pt",
        )
        .Define(
            "Tag_rawFactor_temp", "selPhoton_seedGain[Tag_idx_temp] == 1 ? 1./1.017 : 1"
        )
        .Define("Tag_eta", f"selPhoton_eta[Tag_idx_temp]")
        .Define("Tag_phi", f"selPhoton_phi[Tag_idx_temp]")
        .Define("Tag_mass", "0.0")
        .Define("Tag_label", "2")
        .Define("PhTag_polVec_temp", "ROOT::Math::Polar2DVector(Tag_pt, Tag_phi)")
        .Define("RawPhTag_polVec_temp", "Tag_rawFactor_temp * PhTag_polVec_temp")
    )

    # Probe jet
    rdf = (
        rdf.Define(
            "Jet_indices_temp",
            "findJetIdxs(Jet_pt, Jet_eta, Jet_phi, goodPhoton_eta, goodPhoton_phi, goodPhoton_jetIdx)",
        )
        .Define("Probe_idx_temp", "Jet_indices_temp.first")
        .Filter("Probe_idx_temp >= 0", "Jet found")
        .Filter("Jet_pt[Probe_idx_temp] > 15", "Jet pT > 15")
        .Filter("Jet_vetoed[Probe_idx_temp] == 0", "Probe not vetoed")
        .Filter(
            "fabs(ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)) > 2.7",
            "|dPhi(Probe, Photon)| > 2.7",
        )
        .Filter("Jet_jetId[Probe_idx_temp] >= 4", "Leading jet jet Id >= 4")
        .Define("Activity_idx_temp", "Jet_indices_temp.second")
        .Filter("Probe_idx_temp >= 0", "Jet found")
        .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    rdf = rdf.Define("Activity_denom", "Tag_pt")

    # Redefine T1MET from scratch
    rdf = (
        rdf.Define(
            "metJetIdx_temp", "metJetIdx(Jet_pt, Jet_eta, Jet_phi, Tag_eta, Tag_phi)"
        )
        .Define("metJet_pt_temp", "ROOT::VecOps::Take(Jet_pt, metJetIdx_temp)")
        .Define("metJet_phi_temp", "ROOT::VecOps::Take(Jet_phi, metJetIdx_temp)")
        .Define(
            "metJet_rawFactor_temp", "ROOT::VecOps::Take(Jet_rawFactor, metJetIdx_temp)"
        )
        .Define(
            "metJet_polVec_temp",
            "ROOT::VecOps::Construct<ROOT::Math::Polar2DVector>(metJet_pt_temp, metJet_phi_temp)",
        )
        .Redefine(
            "metJet_polVec_temp",
            "ROOT::VecOps::Sum(metJet_polVec_temp, ROOT::Math::Polar2DVector())",
        )
        .Define(
            "metJet_rawPolVec_temp",
            "ROOT::VecOps::Construct<ROOT::Math::Polar2DVector>((1.0 - metJet_rawFactor_temp) * metJet_pt_temp, metJet_phi_temp)",
        )
        .Redefine(
            "metJet_rawPolVec_temp",
            "ROOT::VecOps::Sum(metJet_rawPolVec_temp, ROOT::Math::Polar2DVector())",
        )
        .Define(
            "PhT1MET_polar_temp",
            "ROOT::Math::Polar2DVector(RawPuppiMET_pt, RawPuppiMET_phi)",
        )
        .Redefine(
            "PhT1MET_polar_temp",
            "PhT1MET_polar_temp - metJet_polVec_temp + metJet_rawPolVec_temp - RawPhTag_polVec_temp + PhTag_polVec_temp",
        )
        .Redefine("T1MET_pt", "float(PhT1MET_polar_temp.R())")
        .Redefine("T1MET_phi", "float(PhT1MET_polar_temp.Phi())")
    )

    for column in jet_columns:
        rdf = rdf.Define("Probe_" + column[4:], f"{column}[Probe_idx_temp]")

    # PhotonJet for cuts
    rdf = rdf.Define(
        "PhotonJet_temp",
        "selPhoton_jetIdx[0] < 0 ? ROOT::Math::Polar2DVector(0, 0) : \
                ROOT::Math::Polar2DVector(Jet_pt[selPhoton_jetIdx[0]], Jet_phi[selPhoton_jetIdx[0]])",
    ).Define(
        "PhotonJet_raw_temp",
        "(1.0 - Jet_rawFactor[selPhoton_jetIdx[0]]) * PhotonJet_temp",
    )

    rdf = (
        rdf.Filter("fabs(1.0-Probe_pt / Tag_pt) < 0.7", "|1-DB| < 0.7; gamjet pass_bal")
        .Filter(
            "fabs(PhT1MET_polar_temp.Dot(PhTag_polVec_temp)/(Tag_pt*Tag_pt)) < 0.7",
            "|1-MPF| < 0.7; gamjet pass_mpf",
        )
        .Filter(
            "selPhoton_jetIdx[0] >= 0 ? \
                        ROOT::VecOps::DeltaR(Jet_eta[selPhoton_jetIdx[0]], Tag_eta, Jet_phi[selPhoton_jetIdx[0]], Tag_phi) < 0.4 : 1",
            "PhotonJet dR < 0.4",
        )
        .Filter(
            "selPhoton_jetIdx[0] >= 0 ? (PhotonJet_raw_temp - RawPhTag_polVec_temp).R() < 0.06*Tag_pt : 1",
            "PhotonJet pt < 0.06*Tag pt",
        )
    )

    return rdf
