#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "./JetMETObjects/interface/FactorizedJetCorrector.h"
#include "./JetMETObjects/interface/JetCorrectorParameters.h"
#include "./JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "./JetMETObjects/interface/FactorizedJetCorrectorCalculator.h"
#include "./JetMETObjects/interface/JetResolutionObject.h"
#include "./JetMETObjects/interface/JetResolution.h"

#include "TString.h"

std::vector<FactorizedJetCorrector*> jet_correctors;
JetCorrectorParameters *L1JetPar;
JetCorrectorParameters* L2RelativeJetPar;
JetCorrectorParameters *L2L3JetPar;

JME::JetResolution *jet_resolution;
JME::JetResolutionScaleFactor *jet_resolution_sf;

// std::vector<JetCorrectorParameters>> allCorrectionParameters;
// JME::JetResolution *jet_resolution(0);
// JetResolutionScaleFactor *jet_resolution_sf(0);

void init_JEC(std::string L1 = "", std::string L2Relative = "", std::string L2L3 = "", unsigned int nThreads = 1) {

    for (unsigned int i = 0; i < nThreads; i++) {
        std::vector<JetCorrectorParameters> correctionParameters(0);
        if (L1 != "") {
            const char *L1_char = L1.c_str();
            L1JetPar = new JetCorrectorParameters(L1_char);
            correctionParameters.push_back(*L1JetPar);
        }
        if (L2Relative != "") {
            const char *L2Relative_char = L2Relative.c_str();
            L2RelativeJetPar = new JetCorrectorParameters(L2Relative_char);
            correctionParameters.push_back(*L2RelativeJetPar);
        }
        if (L2L3 != "") {
            const char *L2L3_char = L2L3.c_str();
            L2L3JetPar = new JetCorrectorParameters(L2L3_char);
            correctionParameters.push_back(*L2L3JetPar);
        }
        if (correctionParameters.size() == 0) {
            std::cout << "JEC files not found" << std::endl;
            return;
        }
        jet_correctors.push_back(new FactorizedJetCorrector(correctionParameters));
    }
}

void init_JER(std::string Res = "", std::string SF = "") {
    // JetCorrectorParameters *resolutionParameters = new JetCorrectorParameters(Res);
    jet_resolution = new JME::JetResolution(Res.c_str());
    if (SF != "") {
        // JetCorrectorParameters *resolution_sf_parameters = new JetCorrectorParameters(SF);
        jet_resolution_sf = new JME::JetResolutionScaleFactor(SF.c_str());
    }
}

ROOT::RVec<float> getJEC(unsigned int threadNumber, ROOT::RVec<float> pt, ROOT::RVec<float> eta, ROOT::RVec<float> area, float rho) {
    ROOT::RVec<float> jec(0);
    for (unsigned int i = 0; i < pt.size(); i++) {
        jet_correctors[threadNumber]->setJetPt(pt[i]);
        jet_correctors[threadNumber]->setJetEta(eta[i]);
        jet_correctors[threadNumber]->setJetA(area[i]);
        jet_correctors[threadNumber]->setRho(rho);
        jec.push_back(jet_correctors[threadNumber]->getCorrection());
    }
    return jec;
}

ROOT::RVec<float> getJER(ROOT::RVec<float> pt, ROOT::RVec<float> eta, ROOT::RVec<float> genJetPt, ROOT::RVec<float> jetArea, ROOT::RVec<float> rho, ROOT::RVec<int> nPU, ROOT::RVec<int> nPV, ROOT::RVec<int> genJetIdx) {
    ROOT::RVec<float> jer(0);
    for (unsigned int i = 0; i < pt.size(); i++) {
        jer.push_back(jet_resolution->getResolution({{JME::Binning::JetPt, pt[i]}, {JME::Binning::JetEta, eta[i]}, {JME::Binning::Rho, rho[i]}}));
    }
    return jer;
}

ROOT::RVec<float> getJER_SF(ROOT::RVec<float> pt, ROOT::RVec<float> eta, ROOT::RVec<float> genJetPt, ROOT::RVec<float> jetArea, ROOT::RVec<float> rho, ROOT::RVec<int> nPU, ROOT::RVec<int> nPV, ROOT::RVec<int> genJetIdx) {
    ROOT::RVec<float> jer_sf(0);
    for (unsigned int i = 0; i < pt.size(); i++) {
        jer_sf.push_back(jet_resolution_sf->getScaleFactor({{JME::Binning::JetPt, pt[i]}, {JME::Binning::JetEta, eta[i]}, {JME::Binning::Rho, rho[i]}, {JME::Binning::NPV, nPV[i]}}));
    }
    return jer_sf;
}
