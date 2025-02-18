import ROOT

def correct_jetId(rdf):
    ROOT.gInterpreter.Declare("""
    #ifndef JETID_FIX
    #define JETID_FIX

    
ROOT::VecOps::RVec<unsigned int> fix_jetId(const ROOT::VecOps::RVec<int>& Jet_jetId, const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_neHEF, const ROOT::VecOps::RVec<float>& Jet_neEmEF, const ROOT::VecOps::RVec<float>& Jet_muEF, const ROOT::VecOps::RVec<float>& Jet_chEmEF) {
    ROOT::VecOps::RVec<unsigned int> Jet_passJetId(Jet_eta.size(), 0);
    
    for (size_t i = 0; i < Jet_eta.size(); ++i) {
        bool passTight = false;
        bool passTightLepVeto = false;
        
        if (abs(Jet_eta[i]) <= 2.7)
            passTight = Jet_jetId[i] & (1 << 1);
        else if (abs(Jet_eta[i]) > 2.7 && abs(Jet_eta[i]) <= 3.0)
            passTight = (Jet_jetId[i] & (1 << 1)) && (Jet_neHEF[i] < 0.99);
        else if (abs(Jet_eta[i]) > 3.0)
            passTight = (Jet_jetId[i] & (1 << 1)) && (Jet_neEmEF[i] < 0.4);
        
        if (passTight) {
            passTightLepVeto = (abs(Jet_eta[i]) > 2.7) ? passTight : (passTight && (Jet_muEF[i] < 0.8) && (Jet_chEmEF[i] < 0.8));
            Jet_passJetId[i] = passTightLepVeto ? 6 : 2;
        }
    }
    return Jet_passJetId;
}

ROOT::VecOps::RVec<unsigned int> fix_jetId(const ROOT::VecOps::RVec<float>& Jet_eta, const ROOT::VecOps::RVec<float>& Jet_neHEF, const ROOT::VecOps::RVec<float>& Jet_neEmEF, const ROOT::VecOps::RVec<int>& Jet_chMultiplicity, const ROOT::VecOps::RVec<int>& Jet_neMultiplicity, const ROOT::VecOps::RVec<float>& Jet_chHEF, const ROOT::VecOps::RVec<float>& Jet_muEF, const ROOT::VecOps::RVec<float>& Jet_chEmEF) {
    ROOT::VecOps::RVec<unsigned int> Jet_passJetId(Jet_eta.size(), 0);
    
    for (size_t i = 0; i < Jet_eta.size(); ++i) {
        bool passTight = false;
        bool passTightLepVeto = false;
        
        if (abs(Jet_eta[i]) <= 2.6)
            passTight = (Jet_neHEF[i] < 0.99) && (Jet_neEmEF[i] < 0.9) && (Jet_chMultiplicity[i] + Jet_neMultiplicity[i] > 1) && (Jet_chHEF[i] > 0.01) && (Jet_chMultiplicity[i] > 0);
        else if (abs(Jet_eta[i]) > 2.6 && abs(Jet_eta[i]) <= 2.7)
            passTight = (Jet_neHEF[i] < 0.90) && (Jet_neEmEF[i] < 0.99);
        else if (abs(Jet_eta[i]) > 2.7 && abs(Jet_eta[i]) <= 3.0)
            passTight = (Jet_neHEF[i] < 0.99);
        else if (abs(Jet_eta[i]) > 3.0)
            passTight = (Jet_neMultiplicity[i] >= 2) && (Jet_neEmEF[i] < 0.4);
        
        if (passTight) {
            passTightLepVeto = (abs(Jet_eta[i]) > 2.7) ? passTight : (passTight && (Jet_muEF[i] < 0.8) && (Jet_chEmEF[i] < 0.8));
            Jet_passJetId[i] = passTightLepVeto ? 6 : 2;
        }
    }
    return Jet_passJetId;
}
#endif
    """)

    jetcols = [str(col) for col in rdf.GetColumnNames() if str(col).startswith("Jet_")]

    if "Jet_chMultiplicity" in jetcols:
        id_input = "Jet_eta,Jet_neHEF,Jet_neEmEF,Jet_chMultiplicity,Jet_neMultiplicity,Jet_chHEF,Jet_muEF,Jet_chEmEF"
    else:
        id_input = "Jet_jetId,Jet_eta,Jet_neHEF,Jet_neEmEF,Jet_muEF,Jet_chEmEF"

    rdf = rdf.Redefine("Jet_jetId", f"fix_jetId({id_input})")

    return rdf

def correct_jets(rdf, cfile, cstack, ccols="Jet_rawPt,Jet_eta,Rho_fixedGridRhoFastjetAll,Jet_area,Jet_phi"):
    ROOT.gInterpreter.Declare(
f"""
#ifndef JET_CORRECTIONS
#define JET_CORRECTIONS

#include <string>
#include "correction.h"

auto cset = correction::CorrectionSet::from_file("{cfile}");
auto cstack = (cset->compound()).at("{cstack}");
#endif
"""
    )

    ROOT.gInterpreter.Declare(
"""
#ifndef GET_CORRECTION
#define GET_CORRECTION
const ROOT::VecOps::RVec<float> get_correction( const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, float rho, const ROOT::VecOps::RVec<float>& area) {
    ROOT::VecOps::RVec<float> correction_factors;

    for (size_t i = 0; i < pt.size(); i++) {
        correction_factors.push_back(cstack->evaluate({area[i], eta[i], pt[i], rho}));
    }

    return correction_factors;
}

const ROOT::VecOps::RVec<float> get_correction( const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, float rho, const ROOT::VecOps::RVec<float>& area, const ROOT::VecOps::RVec<float>& phi) {
    ROOT::VecOps::RVec<float> correction_factors;

    for (size_t i = 0; i < pt.size(); i++) {
        correction_factors.push_back(cstack->evaluate({area[i], eta[i], phi[i], pt[i], rho}));
    }

    return correction_factors;
}
#endif
""")

    rdf = (rdf.Define("Jet_rawPt", "(1.0 - Jet_rawFactor) * Jet_pt")
            .Define("Jet_correctionFactor", f"get_correction({ccols})")
            .Redefine("Jet_pt", "Jet_pt * Jet_correctionFactor")
            .Redefine("Jet_mass", "(1.0 - Jet_rawFactor) * Jet_mass * Jet_correctionFactor")
            .Redefine("Jet_rawFactor", "1.0-1.0/Jet_correctionFactor")
            .Define("UncorrectedJet_temp", "ROOT::VecOps::Construct<ROOT::Math::Polar2DVectorF>(Jet_rawPt, Jet_phi)")
            .Redefine("UncorrectedJet_temp", "ROOT::VecOps::Sum(UncorrectedJet_temp, ROOT::Math::Polar2DVectorF())")
            .Define("CorrectedJet_temp", "ROOT::VecOps::Construct<ROOT::Math::Polar2DVectorF>(Jet_pt, Jet_phi)")
            .Redefine("CorrectedJet_temp", "ROOT::VecOps::Sum(CorrectedJet_temp, ROOT::Math::Polar2DVectorF())")
            .Define("PuppiMET_temp", "ROOT::Math::Polar2DVectorF(PuppiMET_pt, PuppiMET_phi)")
            .Define("T1MET_temp", "PuppiMET_temp + UncorrectedJet_temp - CorrectedJet_temp") # -RC?
            .Redefine("PuppiMET_pt", "float(PuppiMET_temp.R())")
            .Redefine("PuppiMET_phi", "float(PuppiMET_temp.Phi())")
    )

    return rdf

def find_vetojets(rdf, vfile, vset, vcols=["Jet_eta", "Jet_phi"]):
    ROOT.gInterpreter.Declare(
f"""
#ifndef VETO_JETS
#define VETO_JETS

#include <string>
#include "correction.h"

auto vset = correction::CorrectionSet::from_file("{vfile}");
auto veval = vset->at("{vset}");
#endif
"""
    )

    ROOT.gInterpreter.Declare(
"""
#ifndef GET_VETO
#define GET_VETO

ROOT::VecOps::RVec<bool> get_veto( const ROOT::VecOps::RVec<float>& eta, 
                                   const ROOT::VecOps::RVec<float>& phi, 
                                   std::string type="jetvetomap") {
    ROOT::VecOps::RVec<bool> veto_flags;

    for (size_t i = 0; i < eta.size(); i++) {
        if (abs(eta[i]) > 5.131 || abs(phi[i]) > 3.14159) {
            veto_flags.push_back(true);
            continue;
        } 
        float veto = veval->evaluate({type, eta[i], phi[i]});
        veto_flags.push_back(veto > 0.0);
    }

    return veto_flags;
}
#endif
"""
    )

    rdf = (rdf.Define("Jet_vetoed", f"get_veto({','.join(vcols)})"))

    return rdf


def filter_json(rdf, filter_json):
    ROOT.gInterpreter.Declare(
"""
#ifndef JSONFILTER
#define JSONFILTER

#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>
#include <string>

using json = nlohmann::json;

json filter_json;

void init_json(std::string jsonFile) {
    std::cout << "Initializing JSON file" << std::endl;
    std::ifstream f(jsonFile);
    filter_json = json::parse(f);
}

bool isGoodLumi(int run, int lumi) {
   for (auto& lumiRange : filter_json[std::to_string(run)]) {
       if (lumi >= lumiRange[0] && lumi <= lumiRange[1]) {
           return true;
       }
   }

    return false;
}

#endif
"""
    )
    ROOT.init_json(filter_json)
    print("Applying golden JSON cut")
    print(f"JSON file: {filter_json}")
    rdf = (rdf.Filter("isGoodLumi(run, luminosityBlock)", "JSON filter"))
    return rdf

def get_Flags(campaign=None):
    # TODO: Implement campaign-specific flags
    flags = [
            "Flag_goodVertices",
            "Flag_globalSuperTightHalo2016Filter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter",
            "Flag_BadPFMuonDzFilter",
            "Flag_hfNoisyHitsFilter",
            "Flag_eeBadScFilter",
            "Flag_ecalBadCalibFilter"
    ]

    return flags

def sort_jets(rdf, jet_columns):
    # Sort jets by pt
    rdf = rdf.Define("Jet_pt_index", "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(Jet_pt))")
    for col in jet_columns:
        rdf = rdf.Redefine(f"{col}", f"ROOT::VecOps::Take({col}, Jet_pt_index)")
    return rdf