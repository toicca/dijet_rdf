import ROOT

def init_photonjet(rdf, jet_columns, state):
    path = state.module_dir
    path = path / "selections" / "photonjet"
    cpp_path = path / "photonjet.cpp"
    so_path = path / "photonjet_cpp.so"
    h_path = path / "photonjet.h"

    # Compile and load the C++ code
    ROOT.gInterpreter.ProcessLine(f'.L {cpp_path}+')
    ROOT.gSystem.Load(str(so_path))
    ROOT.gInterpreter.Declare(f'#include "{h_path}"')

    # Good photons
    rdf = (rdf.Define("isGoodPhoton", "Photon_cutBased == 3 && \
                      Photon_hoe < 0.02148 && Photon_r9 > 0.94 && Photon_r9 < 1.00")
            .Define("goodPhoton_pt", "Photon_pt[isGoodPhoton]")
            .Define("goodPhoton_eta", "Photon_eta[isGoodPhoton]")
            .Define("goodPhoton_phi", "Photon_phi[isGoodPhoton]")
            .Define("goodPhoton_jetIdx", "Photon_jetIdx[isGoodPhoton]")
            .Filter("goodPhoton_pt.size() > 0", "Good photon found")
    )

    # Trigger selected photons
    rdf = (rdf.Define("hasTrg_temp", "hasTrgObj(goodPhoton_eta, goodPhoton_phi, TrigObj_eta, TrigObj_phi, TrigObj_id)")
            .Define("selPhoton_pt", "goodPhoton_pt[hasTrg_temp]")
            .Define("selPhoton_eta", "goodPhoton_eta[hasTrg_temp]")
            .Define("selPhoton_phi", "goodPhoton_phi[hasTrg_temp]")
            .Define("selPhoton_jetIdx", "goodPhoton_jetIdx[hasTrg_temp]")
            .Filter("selPhoton_pt.size() > 0", "Trigger selected photon found")
            .Filter("selPhoton_pt[0] > 15", "Trigger selected photon pT > 15")
            .Filter("abs(selPhoton_eta[0]) < 1.3", "Trigger selected photon |eta| < 1.3")
            .Define("Tag_idx_temp", "0")
            .Define("Tag_pt", f"selPhoton_pt[Tag_idx_temp]")
            .Define("Tag_rawPt", "Tag_pt")
            .Define("Tag_eta", f"selPhoton_eta[Tag_idx_temp]")
            .Define("Tag_phi", f"selPhoton_phi[Tag_idx_temp]")
            .Define("Tag_mass", "0.0")
            .Define("Tag_label", "2")
    )

    # Probe jet
    rdf = (rdf.Define("Jet_indices_temp", "findJetIdxs(Jet_eta, Jet_phi, goodPhoton_eta, goodPhoton_phi)")
            .Define("Probe_idx_temp", "Jet_indices_temp.first")
            .Filter("Probe_idx_temp >= 0", "Jet found")
            .Filter("Jet_pt[Probe_idx_temp] > 12", "Jet pT > 12")
            .Filter("Jet_vetoed[Probe_idx_temp] == 0", "Probe not vetoed")
            .Filter("abs(ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)) > 2.7", "dPhi(Probe, Photon) > 2.7")
            .Filter("Jet_jetId[Probe_idx_temp] >= 4", "Leading jet jet Id >= 4")
            .Define("Activity_idx_temp", "Jet_indices_temp.second")
            .Filter("Probe_idx_temp >= 0", "Jet found")
            .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    for column in jet_columns:
        rdf = rdf.Define("Probe_"+column[4:], f"{column}[Probe_idx_temp]")

    return rdf
