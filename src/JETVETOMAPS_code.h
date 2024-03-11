#include <iostream>
#include <TFile.h>
#include <TH2D.h>
#include <ROOT/RVec.hxx>  // For ROOT::RVec

// Global TH2D vetomap
TH2D *vetomap = nullptr;

void init_vetomap(TString VetoMapFile) {
    TFile *f = TFile::Open(VetoMapFile);
    if (!f) {
        std::cout << "VetoMap file not found" << std::endl;
        return;
    }
    vetomap = (TH2F*)f->Get("jetvetomap");
    if (!vetomap) {
        std::cout << "VetoMap not found" << std::endl;
        return;
    }
    vetomap->SetDirectory(0);
    f->Close();
}

ROOT::RVec<int> isGoodVeto(ROOT::RVec<float> eta, ROOT::RVec<float> phi) {
    ROOT::RVec<int> passed_veto(0);
    for (unsigned int i = 0; i < eta.size(); i++) {
        int bin = VetoMap->FindBin(eta[i], phi[i]);
        int content = VetoMap->GetBinContent(bin);
        if (content > 0) {
            passed_veto.push_back(0);
        }
        else {
            passed_veto.push_back(1);
        }
    }
    return passed_veto;
}