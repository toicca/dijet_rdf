#!/usr/bin/env python

import ROOT
import glob
import os

def compile_JEC():
    ROOT.gInterpreter.AddIncludePath("./src/JetMETObjects/interface")
    ROOT.gInterpreter.AddIncludePath("./src/JetMETObjects/src")

    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/JetCorrectorParameters.cc+")
    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/SimpleJetCorrector.cc+")
    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/FactorizedJetCorrector.cc+")
    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/FactorizedJetCorrectorCalculator.cc+")

    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc+")
    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/JetCorrectionUncertainty.cc+")

    # For JER
    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/JetResolutionObject.cc++")
    ROOT.gROOT.ProcessLine(".L src/JetMETObjects/src/JetResolution.cc++")

def load_JEC():
    ROOT.gInterpreter.AddIncludePath("./src/JetMETObjects/interface/")
    ROOT.gInterpreter.AddIncludePath("./src/JetMETObjects/src/")
    ROOT.gInterpreter.AddIncludePath("src/JetMETCorrections/Modules/interface")
    ROOT.gInterpreter.AddIncludePath("src/JetMETCorrections/Modules/src")

    ROOT.gInterpreter.ProcessLine('#include "FactorizedJetCorrector.h"')
    ROOT.gInterpreter.ProcessLine('#include "FactorizedJetCorrectorCalculator.h"')
    ROOT.gInterpreter.ProcessLine('#include "SimpleJetCorrector.h"')
    ROOT.gInterpreter.ProcessLine('#include "JetCorrectorParameters.h"')
    ROOT.gInterpreter.ProcessLine('#include "JetCorrectionUncertainty.h"')
    ROOT.gInterpreter.ProcessLine('#include "SimpleJetCorrectionUncertainty.h"')
    ROOT.gInterpreter.ProcessLine('#include "JetResolution.h"')
    ROOT.gInterpreter.ProcessLine('#include "JetResolutionObject.h"')

    
    ROOT.gSystem.Load("src/JetMETObjects/src/JetCorrectorParameters_cc.so")
    ROOT.gSystem.Load("src/JetMETObjects/src/SimpleJetCorrector_cc.so")
    ROOT.gSystem.Load("src/JetMETObjects/src/FactorizedJetCorrector_cc.so")
    ROOT.gSystem.Load("src/JetMETObjects/src/FactorizedJetCorrectorCalculator_cc.so")

    ROOT.gSystem.Load("src/JetMETObjects/src/SimpleJetCorrectionUncertainty_cc.so")
    ROOT.gSystem.Load("src/JetMETObjects/src/JetCorrectionUncertainty_cc.so")

    # For JER
    ROOT.gSystem.Load("./src/JetMETObjects/src/JetResolutionObject_cc.so")
    ROOT.gSystem.Load("./src/JetMETObjects/src/JetResolution_cc.so")

def clean_JEC():
    # Remove .so, .d, .pcm files
    files = glob.glob("src/JetMETObjects/src/*.so")
    files += glob.glob("src/JetMETObjects/src/*.d")
    files += glob.glob("src/JetMETObjects/src/*.pcm")
    for file in files:
        os.remove(file)

if __name__ == "__main__":
    clean_JEC()
    compile_JEC()
    load_JEC()