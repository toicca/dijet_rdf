import ROOT


def init_multijet(rdf, jet_columns, state):

    path = state.module_dir
    path = path / "selections" / "multijet"
    cpp_path = path / "multijet.cpp"
    so_path = path / "multijet_cpp.so"
    h_path = path / "multijet.h"

    # Compile and load the C++ code
    ROOT.gSystem.Load(str(so_path))
    ROOT.gInterpreter.Declare(f'#include "{h_path}"')

    # Multi-jet selection
    rdf = (
        rdf.Filter("nJet > 2", "nJet > 2")
        .Define(
            "RecoilJet_idx_temp",
            "findRecoilJetIdxs(Jet_pt, Jet_eta, Jet_phi, Jet_mass)",
        )
        .Define("nRecoilJet", "int(RecoilJet_idx_temp.size())")
        # .Define("RecoilJet_vetoed", "ROOT::VecOps::Take(Jet_vetoed, RecoilJet_idx_temp)")
        .Filter("RecoilJet_idx_temp.size() >= 2", "At least two recoil jets after veto")
        .Filter("multijetVetoForward(Jet_pt, Jet_eta)", "No jets in |eta| >= 2.5")
        .Filter("multijetVetoNear(Jet_pt, Jet_eta, Jet_phi)", "No jets near lead jet")
        .Filter(
            "Jet_pt[0] > 30 && fabs(Jet_eta[0]) < 1.3 && Jet_jetId[0] >= 4",
            "Leading jet pT > 30, |eta| < 1.3, jetId >= 4",
        )
        .Filter("Jet_vetoed[0] == 0", "Lead jet not vetoed")
        .Filter(
            "Jet_pt[1] > 30 && fabs(Jet_eta[1]) < 2.5 && Jet_jetId[1] >= 4",
            "Second jet pT > 30, |eta| < 2.5, jetId >= 4",
        )
        .Filter("Jet_vetoed[1] == 0", "Second jet not vetoed")
        .Filter(
            "Jet_pt[2] > 30 && fabs(Jet_eta[2]) < 2.5 && Jet_jetId[2] >= 4",
            "Third jet pT > 30, |eta| < 2.5, jetId >= 4",
        )
        .Filter("Jet_vetoed[2] == 0", "Third jet not vetoed")
        .Define("Probe_pt", "Jet_pt[0]")
        .Define("Probe_eta", "Jet_eta[0]")
        .Define("Probe_phi", "Jet_phi[0]")
        .Define("Probe_mass", "Jet_mass[0]")
        .Define("Probe_isFirst", "true")
        .Define("Tag_label", "3")
    )

    # Recoil jets
    rdf = (
        rdf.Define("RecoilJet_pt", "ROOT::VecOps::Take(Jet_pt, RecoilJet_idx_temp)")
        .Define("RecoilJet_eta", "ROOT::VecOps::Take(Jet_eta, RecoilJet_idx_temp)")
        .Define("RecoilJet_phi", "ROOT::VecOps::Take(Jet_phi, RecoilJet_idx_temp)")
        .Define("RecoilJet_mass", "ROOT::VecOps::Take(Jet_mass, RecoilJet_idx_temp)")
        .Define(
            "TagMJ_fourVec_temp",
            "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(RecoilJet_pt, \
                        RecoilJet_eta, RecoilJet_phi, RecoilJet_mass)",
        )
        .Redefine(
            "TagMJ_fourVec_temp",
            "ROOT::VecOps::Sum(TagMJ_fourVec_temp, ROOT::Math::PtEtaPhiMVector())",
        )
    )

    # Tag definitions
    rdf = (
        rdf.Define("Tag_pt", "float(TagMJ_fourVec_temp.Pt())")
        .Define("Tag_eta", "float(TagMJ_fourVec_temp.Eta())")
        .Define("Tag_phi", "float(TagMJ_fourVec_temp.Phi())")
        .Define("Tag_mass", "float(TagMJ_fourVec_temp.M())")
        .Define("Tag_rawPt", "-1.")  # -1 as a place holder
        .Filter("Jet_pt[1] < 0.6*Tag_pt", "Second jet pT < 0.6*recoil pT")
        .Filter("Jet_pt[2] < 0.6*Tag_pt", "Third jet pT < 0.6*recoil pT")
        .Filter(
            "fabs(ROOT::VecOps::DeltaPhi(Tag_phi, Probe_phi)) > 2.84",
            "|dPhi(T,P)| > 2.84",
        )
        .Define("Activity_idx_temp", "-1")  # No activity jet for multijet
    )

    rdf = rdf.Define("Activity_denom", "1.0")

    # Setting up probe columns
    for column in jet_columns:
        if column in ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"]:
            continue
        rdf = rdf.Define("Probe_" + column[4:], f"{column}[0]")

    # Redefine PuppiMET through the lead and recoil
    # rdf = (rdf.Define("MJPuppiMET_temp", "-1.0 * TagMJ_fourVec_temp \
    # - ROOT::Math::PtEtaPhiMVector(Probe_pt, Probe_eta, Probe_phi, Probe_mass)")
    # .Redefine("PuppiMET_pt", "float(MJPuppiMET_temp.Pt())")
    # .Redefine("PuppiMET_phi", "float(MJPuppiMET_temp.Phi())")
    # )

    return rdf
