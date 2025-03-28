import ROOT


def init_zee(rdf, jet_columns, state):
    path = state.module_dir
    path = path / "selections" / "zee"
    cpp_path = path / "zee.cpp"
    so_path = path / "zee_cpp.so"
    h_path = path / "zee.h"

    # Compile and load the C++ code
    ROOT.gSystem.Load(str(so_path))
    ROOT.gInterpreter.Declare(f'#include "{h_path}"')

    # Good electrons selection
    rdf = (
        rdf.Define(
            "electronMask",
            "Electron_pt > 15 && Electron_cutBased >= 3 && Electron_pfRelIso03_all < 1.0",
        )
        .Define("goodElectron_pt", "Electron_pt[electronMask]")
        .Define("goodElectron_eta", "Electron_eta[electronMask]")
        .Define("goodElectron_charge", "Electron_charge[electronMask]")
        .Define("goodElectron_phi", "Electron_phi[electronMask]")
        .Define("goodElectron_mass", "Electron_mass[electronMask]")
        .Filter(
            "goodElectron_pt.size() > 1 && goodElectron_pt.size() < 4",
            "2-3 tight electrons",
        )
    )

    # Trigger selected electrons
    rdf = (
        rdf.Define(
            "trigMask",
            "hasTrgObj(goodElectron_eta, goodElectron_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)",
        )
        .Define("selElectron_pt", "goodElectron_pt[trigMask]")
        .Define("selElectron_eta", "goodElectron_eta[trigMask]")
        .Define("selElectron_charge", "goodElectron_charge[trigMask]")
        .Define("selElectron_phi", "goodElectron_phi[trigMask]")
        .Define("selElectron_mass", "goodElectron_mass[trigMask]")
        .Filter(
            "selElectron_pt.size() > 1 && selElectron_pt.size() < 4",
            "2-3 tight electrons with trigger match",
        )
        .Define(
            "Electron_idx_temp",
            "findElectronIdxs(selElectron_pt, selElectron_eta, selElectron_phi, selElectron_mass, selElectron_charge)",
        )
        .Filter(
            "Electron_idx_temp.first >= 0 && Electron_idx_temp.second >= 0",
            "Two electrons found",
        )
        .Filter(
            "selElectron_pt[Electron_idx_temp.first] > 25 && selElectron_pt[Electron_idx_temp.second] > 15",
            "Leading electron pT > 20, subleading electron pT > 12",
        )
        .Filter(
            "abs(selElectron_eta[Electron_idx_temp.first]) <= 2.4 && abs(selElectron_eta[Electron_idx_temp.second]) <= 2.4",
            "Electron eta <= 2.3",
        )
    )

    # Z boson reconstruction / Tag object
    rdf = (
        rdf.Define(
            "Z_4vec_temp",
            "ROOT::Math::PtEtaPhiMVector(selElectron_pt[Electron_idx_temp.first], \
                        selElectron_eta[Electron_idx_temp.first], selElectron_phi[Electron_idx_temp.first], \
                        selElectron_mass[Electron_idx_temp.first]) + \
                        ROOT::Math::PtEtaPhiMVector(selElectron_pt[Electron_idx_temp.second], \
                        selElectron_eta[Electron_idx_temp.second], selElectron_phi[Electron_idx_temp.second], \
                        selElectron_mass[Electron_idx_temp.second])",
        )
        # Tag definitions
        .Define("Tag_pt", "static_cast<float>(Z_4vec_temp.Pt())")
        .Define("Tag_rawPt", "Tag_pt")
        .Define("Tag_eta", "static_cast<float>(Z_4vec_temp.Eta())")
        .Define("Tag_phi", "static_cast<float>(Z_4vec_temp.Phi())")
        .Define("Tag_mass", "static_cast<float>(Z_4vec_temp.M())")
        .Define("Tag_label", "1")
        .Filter("Tag_pt > 12", "Z pT > 12")
        .Filter("Tag_mass > 71.1876 && Tag_mass < 111.1876", "Z mass window")
    )

    # Probe jet selection
    rdf = (
        rdf.Define(
            "JetElectron_idx_temp",
            "findJetIdxs(Jet_eta, Jet_phi, selElectron_eta, selElectron_phi)",
        )
        .Define("Probe_idx_temp", "JetElectron_idx_temp.first")
        # Probe jet selection
        .Filter("Probe_idx_temp >= 0", "Jet found")
        .Filter("Jet_pt[Probe_idx_temp] > 12", "Leading jet pT > 12")
        .Filter("Jet_jetId[Probe_idx_temp] >= 4", "Leading jet Id >= 4")
        .Filter("Jet_vetoed[Probe_idx_temp] == 0", "Jet not vetoed")
        .Filter(
            "abs(ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)) > 2.7",
            "dPhi(Z,jet) > 2.7",
        )
        .Define(
            "Tag_deltaPhi", "ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)"
        )
        .Define("Activity_idx_temp", "JetElectron_idx_temp.second")
        .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    rdf = rdf.Define("Activity_denom", "Tag_pt")

    for column in jet_columns:
        rdf = rdf.Define("Probe_" + column[4:], f"{column}[Probe_idx_temp]")

    return rdf
