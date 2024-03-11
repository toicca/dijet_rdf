import ROOT
import json
import numpy as np
from typing import List

#====================================================
#
# Functions for files and parsing
#
#====================================================

def findFiles(file : str) -> List[str]:
    with open(file) as f:
        return [line.strip() for line in f.readlines() if line.strip().endswith(".root")]

def readTriggerList(file : str) -> List[str]:
    with open(file) as f:
        return [line.strip() for line in f.readlines()]

#====================================================
#
# Bins for histograms
#
#====================================================
def get_bins() -> dict:
    bins = {}
    
    # These pt and eta bins are some JEC bins
    bins["pt"] = {}
    bins["pt"]["bins"] = np.array((1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133,
                    153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737,
                    790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
                    2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
                    4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000), dtype=float)
    bins["pt"]["n"] = len(bins["pt"]["bins"]) - 1

    bins["eta"] = {}
    bins["eta"]["bins"] = np.array((-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314,
    -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566,
    -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522,
    -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609,
    0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74,
    1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839,
    4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191), dtype=float)
    bins["eta"]["n"] = len(bins["eta"]["bins"]) - 1
    
    bins["phi"] = {}
    bins["phi"]["bins"] = np.linspace(-3.1416, 3.1416, 100, dtype=float)
    bins["phi"]["n"] = len(bins["phi"]["bins"]) - 1
    
    bins["mjj"] = {}
    bins["mjj"]["bins"] = np.linspace(200, 10000, 200, dtype=float)
    bins["mjj"]["n"] = len(bins["mjj"]["bins"]) - 1
    
    bins["deltaEta"] = {}
    bins["deltaEta"]["bins"] = np.linspace(0, 10, 100, dtype=float)
    bins["deltaEta"]["n"] = len(bins["deltaEta"]["bins"]) - 1
    
    bins["deltaR"] = {}
    bins["deltaR"]["bins"] = np.linspace(0, 10, 100, dtype=float)
    bins["deltaR"]["n"] = len(bins["deltaR"]["bins"]) - 1
    
    bins["response"] = {}
    bins["response"]["bins"] = np.linspace(0, 2, 100, dtype=float)
    bins["response"]["n"] = len(bins["response"]["bins"]) - 1
    
    return bins