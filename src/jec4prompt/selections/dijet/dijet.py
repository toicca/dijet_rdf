import ROOT


def init_dijet(rdf, jet_columns, state):
    path = state.module_dir
    path = path / "selections" / "dijet"
    cpp_path = path / "dijet.cpp"
    so_path = path / "dijet_cpp.so"
    h_path = path / "dijet.h"

    # Compile and load the C++ code
    ROOT.gSystem.Load(str(so_path))
    ROOT.gInterpreter.Declare(f'#include "{h_path}"')

    # Tag-probe pair selection
    rdf = (
        rdf.Define(
            "TnP_idx_temp", "findTagProbeIdxs(Jet_eta, Jet_pt, Jet_phi, Jet_jetId)"
        )
        .Filter(
            "TnP_idx_temp.first.first >= 0 && TnP_idx_temp.first.second >= 0",
            "Tag and probe found",
        )
        .Define("Tag_idx_temp", "TnP_idx_temp.first.first")
        .Define("Probe_idx_temp", "TnP_idx_temp.first.second")
        .Filter("Jet_vetoed[Tag_idx_temp] == 0", "Tag not vetoed")
        .Filter("Jet_vetoed[Probe_idx_temp] == 0", "Probe not vetoed")
    )

    # Tag definitions
    rdf = (
        rdf.Define("Tag_pt", "Jet_pt[Tag_idx_temp]")
        .Define("Tag_eta", "Jet_eta[Tag_idx_temp]")
        .Define("Tag_phi", "Jet_phi[Tag_idx_temp]")
        .Define("Tag_mass", "Jet_mass[Tag_idx_temp]")
        .Define("Tag_rawPt", "(1.0 - Jet_rawFactor[Tag_idx_temp]) * Tag_pt")
        .Define("Tag_label", "0")
        .Define("Activity_idx_temp", "TnP_idx_temp.second")
        .Define("Probe_isFirst", "Probe_idx_temp == 0")
        .Filter(
            "fabs(Jet_pt[Probe_idx_temp] - Tag_pt) / (0.5 * (Jet_pt[Probe_idx_temp] + Tag_pt)) < 10",
            "Probe pt not too far from tag pt",
        )  # Basically redundant
    )

    rdf = rdf.Define(
        "Activity_denom", "0.5 * (Jet_pt[Tag_idx_temp] + Jet_pt[Activity_idx_temp])"
    )

    # Probe definitions
    for column in jet_columns:
        rdf = rdf.Define("Probe_" + column[4:], f"{column}[Probe_idx_temp]")

    return rdf
