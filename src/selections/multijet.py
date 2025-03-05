import ROOT

def init_multijet(rdf, jet_columns):

    ROOT.gInterpreter.Declare("""
    #ifndef MULTIJET_IDXS
    #define MULTIJET_IDXS
                              
    ROOT::RVec<int> findRecoilJetIdxs(const ROOT::RVec<float>& Jet_pt, const ROOT::RVec<float>& Jet_eta,
                              const ROOT::RVec<float>& Jet_phi, const ROOT::RVec<float>& Jet_mass, const ROOT::RVec<int>& Jet_jetId) {
        ROOT::RVec<int> idxs;
                              
        for (int i = 1; i < Jet_pt.size(); i++) {
            if (Jet_pt[i] > 30 && abs(Jet_eta[i]) < 2.5 && Jet_jetId[i] >= 4 && abs(ROOT::VecOps::DeltaPhi(Jet_phi[i], Jet_phi[0])) > 1.0) {
                idxs.push_back(i);
            }
        }
                              
        return idxs;
    }
                              
    #endif
    """)

    rdf = (rdf.Filter("nJet > 3", "nJet > 3")
            .Filter("Jet_pt[0] > 30 && abs(Jet_eta[0]) < 2.5 && Jet_jetId[0] >= 4",
                "Leading jet pT > 30 and |eta| < 2.5")
            .Define("RecoilJet_idx_temp", "findRecoilJetIdxs(Jet_pt, Jet_eta, Jet_phi, Jet_mass, Jet_jetId)")
            .Define("RecoilJet_vetoed", "ROOT::VecOps::Take(Jet_vetoed, RecoilJet_idx_temp)")
            .Redefine("RecoilJet_idx_temp", "RecoilJet_idx_temp[RecoilJet_vetoed == 0]")
            .Filter("RecoilJet_idx_temp.size() >= 2", "At least two recoil jets after veto")
            .Define("Probe_pt", "Jet_pt[0]")
            .Define("Probe_eta", "Jet_eta[0]")
            .Define("Probe_phi", "Jet_phi[0]")
            .Define("Probe_mass", "Jet_mass[0]")
            .Define("Probe_isFirst", "true")
            .Filter("Jet_vetoed[0] == 0", "Lead jet not vetoed")
            .Define("Tag_label", "3")
            .Define("RecoilJet_pt", f"ROOT::VecOps::Take(Jet_pt, RecoilJet_idx_temp)")
            .Define("RecoilJet_eta", f"ROOT::VecOps::Take(Jet_eta, RecoilJet_idx_temp)")
            .Define("RecoilJet_phi", f"ROOT::VecOps::Take(Jet_phi, RecoilJet_idx_temp)")
            .Define("RecoilJet_mass", f"ROOT::VecOps::Take(Jet_mass, RecoilJet_idx_temp)")
            .Define("TagMJ_fourVec_temp",
                "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(RecoilJet_pt, \
                        RecoilJet_eta, RecoilJet_phi, RecoilJet_mass)")
            .Redefine("TagMJ_fourVec_temp",
                "ROOT::VecOps::Sum(TagMJ_fourVec_temp, ROOT::Math::PtEtaPhiMVector())")
            .Define("Tag_pt",
                "float(TagMJ_fourVec_temp.Pt())")
            .Define("Tag_eta",
                "float(TagMJ_fourVec_temp.Eta())")
            .Define("Tag_phi",
                "float(TagMJ_fourVec_temp.Phi())")
            .Define("Tag_mass",
                "float(TagMJ_fourVec_temp.M())")
            .Define("Tag_rawPt", "Tag_pt")
            .Define("Activity_idx_temp", "-1") # No activity jet for multijet
    )

    for column in jet_columns:
        if column in ["Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"]:
            continue
        # For multijet change Probe columns to be zero, as probe is not a jet
        rdf = rdf.Define("Probe_"+column[4:], f"{column}[0]")

    return rdf
