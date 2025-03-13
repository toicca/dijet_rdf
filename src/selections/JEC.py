jet_columns = [
    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId",
    "Jet_area", "Jet_nConstituents", "Jet_nElectrons", "Jet_nMuons",
    "Jet_chEmEF", "Jet_chHEF",
    "Jet_neEmEF", "Jet_neHEF",
    "Jet_hfEmEF", "Jet_hfHEF",
    "Jet_muEF",
    "Jet_neMultiplicity", "Jet_chMultiplicity",
    "Jet_rawFactor",
    "Jet_btagPNetQvG"
]

def _init_TnP(rdf, state):
    """
    Initialize the tag and probe variables for the analysis
    """
    channel = state.args.channel
    logger = state.logger   
    if channel == "dijet":
        from selections.dijet.dijet import init_dijet as init_selection
    elif channel == "zmm":
        from selections.zmm.zmm import init_zmm as init_selection
    elif channel == "zee":
        from selections.zee.zee import init_zee as init_selection
    elif channel == "photonjet":
        from selections.photonjet.photonjet import init_photonjet as init_selection
    elif channel == "multijet":
        from selections.multijet.multijet import init_multijet as init_selection
    else:
        logger.info("NOTICE: Running on empty selection")
        from src.selections.empty.empty import init_empty as init_selection

    rdf = init_selection(rdf, jet_columns, state)

    probe_cols = [str(col) for col in rdf.GetColumnNames() if str(col).startswith("Probe_")]
    if "Probe_rawPt" not in probe_cols:
        rdf = rdf.Define("Probe_rawPt", "(1.0 - Probe_rawFactor) * Probe_pt")
    rdf = rdf.Filter("Activity_idx_temp >= 0 ? (Jet_pt[Activity_idx_temp] / ((Probe_pt + Tag_pt)*0.5)) < 1.0 : 1", "Activity jet pT fraction < 1.0")

    # Label non-flat branches as _temp to drop them later
    rdf = (rdf.Define("Tag_fourVec_temp", "ROOT::Math::PtEtaPhiMVector(Tag_pt, Tag_eta, Tag_phi, Tag_mass)")
            .Define("Probe_fourVec_temp", "ROOT::Math::PtEtaPhiMVector(Probe_pt, Probe_eta, \
                    Probe_phi, Probe_mass)")
            .Define("Tag_polarVec_temp", "ROOT::Math::Polar2DVector(Tag_pt, Tag_phi)")
            .Define("Tag_raw_polarVec_temp", "ROOT::Math::Polar2DVector(Tag_rawPt, Tag_phi)")
            .Define("Probe_polarVec_temp", "ROOT::Math::Polar2DVector(Probe_pt, Probe_phi)")
            .Define("Probe_raw_polarVec_temp", "ROOT::Math::Polar2DVector(Probe_rawPt, Probe_phi)")
            .Define("PuppiMET_polarVec_temp", "ROOT::Math::Polar2DVector(PuppiMET_pt, \
                    PuppiMET_phi)")
            .Define("RawPuppiMET_polarVec_temp", "ROOT::Math::Polar2DVector(RawPuppiMET_pt, \
                    RawPuppiMET_phi)")
            .Define("T1MET_polarVec_temp", "ROOT::Math::Polar2DVector(T1MET_pt, T1MET_phi)")
    )

    # Activity vector for HDM
    # For dijet and multijet this is the sum of all the jets minus the tag and probe
    # For zjet and egamma this is the sum of all the jets minus the probe
    rdf = (rdf.Define("JetActivity_fourVec_temp", 
                    "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Jet_pt, Jet_eta, \
                            Jet_phi, Jet_mass)")
            .Redefine("JetActivity_fourVec_temp", "ROOT::VecOps::Sum(JetActivity_fourVec_temp, \
                    ROOT::Math::PtEtaPhiMVector())")
    )
    if channel == "dijet" or channel == "multijet":
        rdf = (rdf.Redefine("JetActivity_fourVec_temp",
                            "JetActivity_fourVec_temp - Tag_fourVec_temp - Probe_fourVec_temp"))
    elif channel == "zjet" or channel == "egamma":
        rdf = (rdf.Redefine("JetActivity_fourVec_temp",
                            "JetActivity_fourVec_temp - Probe_fourVec_temp"))
    rdf = (rdf.Define("JetActivity_pt", "float(JetActivity_fourVec_temp.Pt())")
            .Define("JetActivity_eta", "float(JetActivity_fourVec_temp.Eta())")
            .Define("JetActivity_phi", "float(JetActivity_fourVec_temp.Phi())")
            .Define("JetActivity_mass", "float(JetActivity_fourVec_temp.M())")
            .Define("JetActivity_polarVec_temp",
                "ROOT::Math::Polar2DVector(JetActivity_pt, JetActivity_phi)")
    )

    rdf = (rdf.Define("TnP_pt_ave", "(Probe_pt + Tag_pt) * 0.5")
            .Define("TnP_deltaPhi", "ROOT::VecOps::DeltaPhi(Tag_phi, Probe_phi)")
            .Define("TnP_deltaR", "ROOT::VecOps::DeltaR(Tag_eta, Probe_eta, Tag_phi, Probe_phi)")
    )
    
    return rdf


def _def_JEC(rdf):
    rdf = (rdf.Define("DB_direct",
                "-1.0 * Tag_polarVec_temp.Dot(Probe_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("DB_raw_direct",
                "-1.0 * Tag_raw_polarVec_temp.Dot(Probe_raw_polarVec_temp) / (Tag_rawPt * Tag_rawPt)")
            .Define("DB_ratio", "Probe_pt / Tag_pt")
            .Define("DB_raw_ratio", "Probe_rawPt / Tag_rawPt")
            .Define("DBT1", "1.0 + T1MET_polarVec_temp.Dot(Tag_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("MPF_tag",
                "1.0 + PuppiMET_polarVec_temp.Dot(Tag_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("MPF_raw_tag",
                "1.0 + RawPuppiMET_polarVec_temp.Dot(Tag_raw_polarVec_temp) / (Tag_rawPt * Tag_rawPt)")
            .Define("MPF_probe",
                "1.0 + PuppiMET_polarVec_temp.Dot(Probe_polarVec_temp) / (Probe_pt * Probe_pt)")
            .Define("MPF_raw_probe",
                "1.0 + RawPuppiMET_polarVec_temp.Dot(Probe_raw_polarVec_temp) / (Probe_rawPt * Probe_rawPt)")
            .Define("R_un_reco_tag_temp",
                "JetActivity_polarVec_temp.Dot(Tag_polarVec_temp) / (Tag_pt * Tag_pt)")
            .Define("R_un_gen_tag_temp", "1.0")
            .Define("R_un_reco_probe_temp",
                "JetActivity_polarVec_temp.Dot(Probe_polarVec_temp) / (Probe_pt * Probe_pt)")
            .Define("R_un_gen_probe_temp", "1.0")
            .Define("HDM_tag",
                "(DB_direct + MPF_tag - 1.0 + R_un_reco_tag_temp - R_un_gen_tag_temp) / \
                        (cos(ROOT::VecOps::DeltaPhi(Tag_phi, Probe_phi)))")
            .Define("HDM_probe",
                "(DB_direct + MPF_probe - 1.0 + R_un_reco_probe_temp - R_un_gen_probe_temp) / \
                        (cos(ROOT::VecOps::DeltaPhi(Tag_phi, Probe_phi)))")
           )


    # Energy Fraction balance
    rdf = (rdf.Define("EFB_chEmEF", "(Probe_rawPt * Probe_chEmEF) / Tag_pt")
        .Define("EFB_chHEF", "(Probe_rawPt * Probe_chHEF) / Tag_pt")
        .Define("EFB_CHF", "(Probe_rawPt * (Probe_chHEF + Probe_chEmEF)) / Tag_pt")
        .Define("EFB_hfEmEF", "(Probe_rawPt * Probe_hfEmEF) / Tag_pt")
        .Define("EFB_hfHEF", "(Probe_rawPt * Probe_hfHEF) / Tag_pt")
        .Define("EFB_NHF", "(Probe_rawPt * (Probe_hfHEF + Probe_hfEmEF)) / Tag_pt")
        .Define("EFB_muEF", "(Probe_rawPt * Probe_muEF) / Tag_pt")
        .Define("EFB_neEmEF", "(Probe_rawPt * Probe_neEmEF) / Tag_pt")
        .Define("EFB_neHEF", "(Probe_rawPt * Probe_neHEF) / Tag_pt")
        .Define("EFB_NEF", "(Probe_rawPt * (Probe_neHEF + Probe_neEmEF)) / Tag_pt")
    )

    return rdf

def _check_JEC(rdf):
    # Define energy fraction variables in case not in NanoAOD (added in V14)
    ef_cols = ["Jet_chEmEF", "Jet_chHEF", "Jet_hfEmEF", "Jet_muEF", "Jet_neEmEF", "Jet_neHEF"]
    rdf_cols = [str(col) for col in rdf.GetColumnNames() if str(col).startswith("Jet_")]

    for ef in ef_cols:
        if ef not in rdf_cols:
            rdf = rdf.Define(ef, "-1.0")


def run_JEC(rdf, state):
    logger = state.logger
    _check_JEC(rdf)
    # Initialize the JEC variables
    logger.info("Initializing TnP variables")
    rdf = _init_TnP(rdf, state)
    logger.info("Initializing JEC variables")
    rdf = _def_JEC(rdf)

    return rdf
