def init_empty(rdf, jet_columns, state):
    rdf = (
        rdf.Define("Tag_pt", "-1.0")
        .Define("Tag_eta", "0.0")
        .Define("Tag_phi", "0.0")
        .Define("Tag_mass", "0.0")
        .Define("Tag_rawPt", "-1.0")
        .Define("Tag_label", "-1")
        .Define("Probe_pt", "-1.0")
        .Define("Probe_eta", "0.0")
        .Define("Probe_phi", "0.0")
        .Define("Probe_mass", "0.0")
        .Define("Probe_rawPt", "-1.0")
        .Define("Activity_idx_temp", "-1")
    )

    for column in jet_columns:
        if column in ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"]:
            continue
        rdf = rdf.Define("Probe_" + column[4:], "0.0")

    return rdf
