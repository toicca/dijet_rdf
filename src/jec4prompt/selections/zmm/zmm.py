import ROOT


def init_zmm(rdf, jet_columns, state):
    ROOT.gInterpreter.Declare("""
        #ifndef ZMM_SEL
        #define ZMM_SEL

        #include "ROOT/RVec.hxx"
        #include <cmath>
        #include "Math/Vector4D.h"

        std::pair<int, int> findMuonIdxs(const ROOT::RVec<float>& Muon_pt,
                                         const ROOT::RVec<float>& Muon_eta,
                                         const ROOT::RVec<float>& Muon_phi,
                                         const ROOT::RVec<float>& Muon_mass,
                                         const ROOT::RVec<float>& Muon_charge) {
            int idx1 = -1;
            int idx2 = -1;

            float Zmass = 91.1876;
            float mtemp = 0.;

            for (int i = 0; i < Muon_pt.size(); i++) {
                ROOT::Math::PtEtaPhiMVector mu1(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
                for (int j = i+1; j < Muon_pt.size(); j++) {
                    if (Muon_charge[i] == Muon_charge[j]) continue;
                    ROOT::Math::PtEtaPhiMVector mu2(Muon_pt[j], Muon_eta[j], \
                            Muon_phi[j], Muon_mass[j]);
                    ROOT::Math::PtEtaPhiMVector Z = mu1 + mu2;

                    if (abs(mtemp - Zmass) > abs(Z.M() - Zmass)) {
                        mtemp = Z.M();
                        idx1 = i;
                        idx2 = j;
                    }
                }
            }

            if (idx1 == -1 || idx2 == -1) return std::make_pair(-1, -1);

            // Sort muons by pt
            if (Muon_pt[idx1] < Muon_pt[idx2]) {
                int temp = idx1;
                idx1 = idx2;
                idx2 = temp;
            }

            return std::make_pair(idx1, idx2);
        }

        ROOT::RVec<bool> hasTrgObj(const ROOT::RVec<float>& Muon_eta,
                                   const ROOT::RVec<float>& Muon_phi,
                                   const ROOT::RVec<float>& trg_eta,
                                   const ROOT::RVec<float>& trg_phi,
                                   const ROOT::RVec<int>& trg_filterBits,
                                   const ROOT::RVec<int>& trg_id) {
            ROOT::RVec<bool> trgObj(Muon_eta.size(), false);
            for (int i = 0; i < Muon_eta.size(); i++) {
                for (int j = 0; j < trg_eta.size(); j++) {
                    if (ROOT::VecOps::DeltaR(Muon_eta[i], trg_eta[j], \
                            Muon_phi[i], trg_phi[j]) < 0.3) {
                        if (trg_id[j] == 13 && (trg_filterBits[j] & (1 << 3))) {
                            trgObj[i] = true;
                            break;
                        }
                    }
                }
            }
            return trgObj;
        }

        std::pair<int, int> findJetIdxs(const ROOT::RVec<float>& Jet_pt,
                                        const ROOT::RVec<float>& Jet_eta,
                                        const ROOT::RVec<float>& Jet_phi,
                                        const ROOT::RVec<float>& Muon_eta,
                                        const ROOT::RVec<float>& Muon_phi) {
            int idx1 = -1;
            int idx2 = -1;
            for (int i = 0; i < Jet_eta.size(); i++) {
                bool badJet = false;
                for (int j = 0; j < Muon_eta.size(); j++) {
                    if (ROOT::VecOps::DeltaR(Jet_eta[i], Muon_eta[j], \
                            Jet_phi[i], Muon_phi[j]) < 0.3) {
                        badJet = true;
                        break;
                    }
                }
                if (!badJet) {
                    if (idx1 == -1) {
                        idx1 = i;
                    } else if (idx2 == -1 && Jet_pt[i] > 15) {
                        idx2 = i;
                        break;
                    }
                }
            }
            return std::make_pair(idx1, idx2);
        }

        ROOT::RVec<bool> passMuMET(const ROOT::RVec<float>& Muon_eta,
                                const ROOT::RVec<float>& Muon_phi,
                                const ROOT::RVec<float>& Jet_pt,
                                const ROOT::RVec<float>& Jet_eta,
                                const ROOT::RVec<float>& Jet_phi,
                                const ROOT::RVec<float>& Jet_rawFactor,
                                const ROOT::RVec<float>& Jet_muonSubtrFactor) {
            ROOT::RVec<bool> passMuMET(Jet_eta.size(), true);

            // dR check
            for (int i = 0; i < Jet_eta.size(); i++) {
                for (int j = 0; j < Muon_eta.size(); j++) {
                    if (ROOT::VecOps::DeltaR(Jet_eta[i], Muon_eta[j], \
                            Jet_phi[i], Muon_phi[j]) < 0.3) {
                        passMuMET[i] = false;
                        break;
                    }
                }
            }

            // pT check
            // Require that the corrected pT of the jet is greater than 15 GeV after muon subtraction
            // TODO: Should the muon be subtracted from the probe jet?
            for (int i = 0; i < Jet_eta.size(); i++) {
                if (!passMuMET[i]) continue;

                float rpt = (1.0 - Jet_rawFactor[i]) * Jet_pt[i];
                float musubtrpt = rpt * (1.0 - Jet_muonSubtrFactor[i]);
                float mupt = rpt - musubtrpt;

                if (Jet_pt[i] - mupt < 15) {
                    passMuMET[i] = false;
                }
            }


            return passMuMET;
        }
        #endif
    """)

    # Good muon selection
    rdf = (
        rdf.Define(
            "muonMask",
            "Muon_pt > 8 && Muon_pfIsoId >= 4 && Muon_pfRelIso04_all < 0.15 && Muon_tightId",
        )
        .Define("goodMuon_pt", "Muon_pt[muonMask]")
        .Define("goodMuon_eta", "Muon_eta[muonMask]")
        .Define("goodMuon_charge", "Muon_charge[muonMask]")
        .Define("goodMuon_phi", "Muon_phi[muonMask]")
        .Define("goodMuon_mass", "Muon_mass[muonMask]")
    )

    # Trigger selected muons
    rdf = (
        rdf.Define(
            "trigMask",
            "hasTrgObj(goodMuon_eta, goodMuon_phi, TrigObj_eta, TrigObj_phi, TrigObj_filterBits, TrigObj_id)"
        )
        .Define("selMuon_pt", "goodMuon_pt[trigMask]")
        .Define("selMuon_eta", "goodMuon_eta[trigMask]")
        .Define("selMuon_charge", "goodMuon_charge[trigMask]")
        .Define("selMuon_phi", "goodMuon_phi[trigMask]")
        .Define("selMuon_mass", "goodMuon_mass[trigMask]")
        .Filter(
            "selMuon_pt.size() > 1 && selMuon_pt.size() < 4",
            "2-3 tight muons with trigger match"
        )
        .Define(
            "Muon_idx_temp",
            "findMuonIdxs(selMuon_pt, selMuon_eta, selMuon_phi, selMuon_mass, selMuon_charge)"
        )
        .Filter(
            "Muon_idx_temp.first >= 0 && Muon_idx_temp.second >= 0", "Two muons found"
        )
        .Filter(
            "selMuon_pt[Muon_idx_temp.first] > 20 && selMuon_pt[Muon_idx_temp.second] > 10",
            "Leading muon pT > 20, subleading muon pT > 10"
        )
        .Filter(
            "fabs(selMuon_eta[Muon_idx_temp.first]) <= 2.3 && fabs(selMuon_eta[Muon_idx_temp.second]) <= 2.3",
            "|Muon eta| <= 2.3"
        )
    )

    # Z boson reconstruction / Tag object
    rdf = (
        rdf.Define(
            "Z_4vec_temp",
            "ROOT::Math::PtEtaPhiMVector(selMuon_pt[Muon_idx_temp.first], \
                        selMuon_eta[Muon_idx_temp.first], selMuon_phi[Muon_idx_temp.first], \
                        selMuon_mass[Muon_idx_temp.first]) + \
                        ROOT::Math::PtEtaPhiMVector(selMuon_pt[Muon_idx_temp.second], \
                        selMuon_eta[Muon_idx_temp.second], selMuon_phi[Muon_idx_temp.second], \
                        selMuon_mass[Muon_idx_temp.second])"
        )
        .Define("Tag_pt", "static_cast<float>(Z_4vec_temp.Pt())")
        .Define("Tag_rawPt", "Tag_pt")
        .Define("Tag_eta", "static_cast<float>(Z_4vec_temp.Eta())")
        .Define("Tag_phi", "static_cast<float>(Z_4vec_temp.Phi())")
        .Define("Tag_mass", "static_cast<float>(Z_4vec_temp.M())")
        .Define("Tag_label", "1")
        .Filter("Tag_pt > 12",
                "Z pT > 12")
        .Filter("Tag_mass > 71.1876 && Tag_mass < 111.1876",
                "Z mass window, +-20 GeV")
    )

    # Recalculate MET
    # Note: only use jets with pT > 15 GeV after muon subtraction and no deltaR < 0.3 to a selected muon
    rdf = (
        rdf.Define("passMuDr", "passMuMET(selMuon_eta, selMuon_phi, Jet_pt, Jet_eta, Jet_phi, Jet_rawFactor, Jet_muonSubtrFactor)")
        .Define("metJets", "passMuDr && Jet_pt > 15")
        .Define("metJet_pt_temp", "Jet_pt[metJets]")
        .Define("metJet_phi_temp", "Jet_phi[metJets]")
        .Define("metJet_rawFactor_temp", "Jet_rawFactor[metJets]")
        .Define(
            "metJet_polVec_temp",
            "ROOT::VecOps::Construct<ROOT::Math::Polar2DVector>(metJet_pt_temp, metJet_phi_temp)"
        )
        .Redefine(
            "metJet_polVec_temp",
            "ROOT::VecOps::Sum(metJet_polVec_temp, ROOT::Math::Polar2DVector())"
        )
        .Define(
            "metJet_rawPolVec_temp",
            "ROOT::VecOps::Construct<ROOT::Math::Polar2DVector>((1.0 - metJet_rawFactor_temp) * metJet_pt_temp, metJet_phi_temp)"
        )
        .Redefine(
            "metJet_rawPolVec_temp",
            "ROOT::VecOps::Sum(metJet_rawPolVec_temp, ROOT::Math::Polar2DVector())"
        )
        .Define(
            "ZmmT1MET_polar_temp",
            "ROOT::Math::Polar2DVector(RawPuppiMET_pt, RawPuppiMET_phi)"
        )
        .Redefine(
            "ZmmT1MET_polar_temp",
            "ZmmT1MET_polar_temp + metJet_rawPolVec_temp - metJet_polVec_temp"
        )
        .Redefine("T1MET_pt", "float(ZmmT1MET_polar_temp.R())")
        .Redefine("T1MET_phi", "float(ZmmT1MET_polar_temp.Phi())")
    )

    # Probe jet selection
    rdf = (
        rdf.Define(
            "JetMuon_idx_temp",
            "findJetIdxs(Jet_pt, Jet_eta, Jet_phi, selMuon_eta, selMuon_phi)" # Could be goodMuons
        )
        .Define("Probe_idx_temp", "JetMuon_idx_temp.first")
        .Filter("Probe_idx_temp >= 0",
                "Jet found")
        .Filter("Jet_pt[Probe_idx_temp] > 12",
                "Leading jet pT > 12")
        # .Filter("(1.0-Jet_rawFactor[Probe_idx_temp])*Jet_pt[Probe_idx_temp] > 12",
        # "Leading jet raw pT > 12")
        .Filter("Jet_jetId[Probe_idx_temp] >= 4",
                "Leading jet Id >= 4")
        .Filter("Jet_vetoed[Probe_idx_temp] == 0",
                "Jet not vetoed")
        .Filter(
            "fabs(ROOT::VecOps::DeltaPhi(Jet_phi[Probe_idx_temp], Tag_phi)) > 2.7",
            "|dPhi(Z,jet)| > 2.7"
        )
        .Define("Activity_idx_temp", "JetMuon_idx_temp.second")
        .Define("Probe_isFirst", "Probe_idx_temp == 0")
    )

    rdf = rdf.Define("Activity_denom", "Tag_pt")

    for column in jet_columns:
        rdf = rdf.Define("Probe_" + column[4:], f"{column}[Probe_idx_temp]")

    return rdf
