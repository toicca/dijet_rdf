import ROOT

def correct_jets(rdf, cfile, cstack, ccols=["Jet_rawPt", "Jet_eta", "Rho_fixedGridRhoFastjetAll", "Jet_area", "Jet_phi"]):
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
            .Define("Jet_correctionFactor", f"get_correction({','.join(ccols)})")
            .Redefine("Jet_pt", "Jet_pt * Jet_correctionFactor")
            .Redefine("Jet_mass", "(1.0 - Jet_rawFactor) * Jet_mass * Jet_correctionFactor")
            .Redefine("Jet_rawFactor", "1.0-Jet_correctionFactor")
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

const ROOT::VecOps::RVec<int> get_veto( const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, std::string type="jetvetomap") {
    ROOT::VecOps::RVec<bool> veto_flags;

    for (size_t i = 0; i < eta.size(); i++) {
        float veto = veval->evaluate({type, eta[i], phi[i]});
        if (veto > 0.0) {
            veto_flags.push_back(true);
        } else {
            veto_flags.push_back(false);
        }
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
    golden_json = json::parse(f);
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