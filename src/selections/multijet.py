import ROOT

def init_multijet(rdf, jet_columns):
    # Change Tag <-> Probe for multijet, since low pt jets better calibrated?
    recoil_filter = "abs(RecoilJet_eta)<2.5 && RecoilJet_pt>30"
    rdf = (rdf.Filter("nJet > 2", "nJet > 2")
            .Filter("Jet_pt[0] > 30",
                "Leading jet pT > 30")
            .Define("Probe_pt", "Jet_pt[0]")
            .Define("Probe_eta", "Jet_eta[0]")
            .Define("Probe_phi", "Jet_phi[0]")
            .Define("Probe_mass", "Jet_mass[0]")
            .Define("Tag_label", "3")
            .Define("RecoilJet_ids",
                "ROOT::VecOps::Drop(ROOT::VecOps::Enumerate(Jet_pt), \
                        ROOT::VecOps::RVec<int>{0})")
            .Define("RecoilJet_pt", f"ROOT::VecOps::Take(Jet_pt, RecoilJet_ids)")
            .Define("RecoilJet_eta", f"ROOT::VecOps::Take(Jet_eta, RecoilJet_ids)")
            .Define("RecoilJet_phi", f"ROOT::VecOps::Take(Jet_phi, RecoilJet_ids)")
            .Define("RecoilJet_mass", f"ROOT::VecOps::Take(Jet_mass, RecoilJet_ids)")
            .Define("Tag_pt", f"RecoilJet_pt[{recoil_filter}]")
            .Define("Tag_rawPt", "0.0") # TODO
            .Define("Tag_eta", f"RecoilJet_eta[{recoil_filter}]")
            .Define("Tag_phi", f"RecoilJet_phi[{recoil_filter}]")
            .Define("Tag_mass", f"RecoilJet_mass[{recoil_filter}]")
            .Define("Tag_ids", f"RecoilJet_ids[{recoil_filter}]")
            .Define("TagMJ_fourVec_temp",
                "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Tag_pt, \
                        Tag_eta, Tag_phi, Tag_mass)")
            .Redefine("TagMJ_fourVec_temp",
                "ROOT::VecOps::Sum(TagMJ_fourVec_temp, ROOT::Math::PtEtaPhiMVector())")
            .Redefine("Tag_pt",
                "float(TagMJ_fourVec_temp.Pt())")
            .Redefine("Tag_eta",
                "float(TagMJ_fourVec_temp.Eta())")
            .Redefine("Tag_phi",
                "float(TagMJ_fourVec_temp.Phi())")
            .Redefine("Tag_mass",
                "float(TagMJ_fourVec_temp.M())")
            .Define("Activity_idx_temp", "-1") # No activity jet for multijet
    )

    for column in jet_columns:
        if column in ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"]:
            continue
        # For multijet change Probe columns to be zero, as probe is not a jet
        rdf = rdf.Define("Probe_"+column[4:], "0.0")

    return rdf
