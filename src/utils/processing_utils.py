import configparser
import json
import numpy as np
import subprocess
from typing import Dict, List

def find_site(file):
    try:
        result = subprocess.run(
            [
                '/cvmfs/cms.cern.ch/common/dasgoclient',
                '--json',
                f'--query=site file={file}'
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )

        # Parse the JSON output
        data = json.loads(result.stdout)
        site_json = data[0]['site'][0] # Random indices

        # Check available sites
        # Don't include test sites or sites outside Europe
        available_sites = [key for key in site_json['states'] if site_json['states'][key] == 'AVAILABLE' and '_Test' not in key and '_Tape' not in key and '_US_' not in key and '_RU_' not in key and '_BR_' not in key and '_KR_' not in key and '_CN_' not in key]

        final_paths = {}
        for site in available_sites:
            final_paths[site] = site_json['rses'][site][0]

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while processing file '{file}': {e.stderr}")
        return None

    except json.JSONDecodeError as e:
        print(f"Failed to parse JSON output for file '{file}': {e}")
        return None

    return final_paths

def find_luminosity():
    pass

def get_bins(fill_range : tuple = (376370, 380100), isMC : bool = False) -> dict:
    bins = {}
    
    # These pt and eta bins are some JEC bins
    bins["pt"] = {}
    bins["pt"]["bins"] = np.array((1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133,
                    153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468, 507, 548, 592, 638, 686, 737,
                    790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,
                    2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
                    4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000), dtype=float)
    bins["pt"]["n"] = len(bins["pt"]["bins"]) - 1

    # Double_t zptbins[] = {12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500}
    bins["zpt"] = {}
    bins["zpt"]["bins"] = np.array((12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500), dtype=float)
    bins["zpt"]["n"] = len(bins["zpt"]["bins"]) - 1

    bins["zeta"] = {}
    bins["zeta"]["bins"] = np.linspace(-5,5,100, dtype=float)
    bins["zeta"]["n"] = len(bins["zeta"]["bins"]) - 1

    bins["eta"] = {}
    bins["eta"]["bins"] = np.array((-5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314,
    -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653, -1.566,
    -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522,
    -0.435, -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609,
    0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.74,
    1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839,
    4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191), dtype=float)
    bins["eta"]["n"] = len(bins["eta"]["bins"]) - 1

    bins["jet1pt"] = {}
    bins["jet1pt"]["bins"] = np.linspace(0, 1750, 1750, dtype=float)
    bins["jet1pt"]["n"] = len(bins["jet1pt"]["bins"]) - 1

    bins["jet1eta"] = {}
    bins["jet1eta"]["bins"] = np.linspace(-6, 6, 120, dtype=float)
    bins["jet1eta"]["n"] = len(bins["jet1eta"]["bins"]) - 1
    
    bins["phi"] = {}
    bins["phi"]["bins"] = np.linspace(-np.pi, np.pi, 73, dtype=float)
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
    
    bins["asymmetry"] = {}
    bins["asymmetry"]["bins"] = np.linspace(-1, 1, 100, dtype=float)
    bins["asymmetry"]["n"] = len(bins["asymmetry"]["bins"]) - 1
    
    bins["runs"] = {}
    if not isMC:
        # bins["runs"]["bins"] = np.linspace(355065, 391370, int((391370-355065) / 10000), dtype=float)
        bins["runs"]["bins"] = np.linspace(fill_range[0], fill_range[1], int((fill_range[1] - fill_range[0])), dtype=float)
        bins["runs"]["n"] = len(bins["runs"]["bins"]) - 1
    else:
        bins["runs"]["bins"] = np.linspace(0, 2, 2, dtype=float)
        bins["runs"]["n"] = len(bins["runs"]["bins"]) - 1
    
    bins["bx"] = {}
    bins["bx"]["bins"] = np.linspace(0, 3564, 3564, dtype=float)
    bins["bx"]["n"] = len(bins["bx"]["bins"]) - 1
    
    bins["lumi"] = {}
    bins["lumi"]["bins"] = np.linspace(0, 1600, 1600, dtype=float)
    bins["lumi"]["n"] = len(bins["lumi"]["bins"]) - 1

    bins["tagger"] = {}
    bins["tagger"]["bins"] = np.linspace(0, 1, 100, dtype=float)
    bins["tagger"]["n"] = len(bins["tagger"]["bins"]) - 1
    
    bins["jec_pt"] = {}
    bins["jec_pt"]["bins"] = np.array((1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,
    97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,
    507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248,
    1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500,
    2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832, 4037, 4252, 4477, 4713,
    4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000), dtype=float)
    bins["jec_pt"]["n"] = len(bins["jec_pt"]["bins"]) - 1

    bins["nestor_tnp_pt"] = {}
    bins["nestor_tnp_pt"]["bins"] = np.array((10, 15, 21, 28, 37, 49,
    59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575,
    638, 737, 846, 967, 1101, 1248,
    1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103, 3450, 4037, 5220), dtype=float)
    bins["nestor_tnp_pt"]["n"] = len(bins["nestor_tnp_pt"]["bins"]) - 1

    bins["nestor_tnp_eta"] = {}
    bins["nestor_tnp_eta"]["bins"] = np.array((-5.191,
    -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489,
    -3.314, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -2.043,
    -1.93, -1.83, -1.74, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218,
    -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435,
    -0.348, -0.261, -0.174, -0.087, 0, 0.087, 0.174, 0.261, 0.348, 0.435,
    0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305,
    1.392, 1.479, 1.566, 1.653, 1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5,
    2.65, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191,
    4.363, 4.538, 4.716, 4.889, 5.191), dtype=float)
    bins["nestor_tnp_eta"]["n"] = len(bins["nestor_tnp_eta"]["bins"]) - 1

    bins["beta"] = {}
    bins["beta"]["bins"] = np.array((-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.650,
                 -2.500, -2.322, -2.172, -1.930, -1.653, -1.479, -1.305,
                 -1.044, -0.783, -0.522, -0.261, 0.000, 0.261, 0.522, 0.783,
                 1.044, 1.305, 1.479, 1.653, 1.930, 2.172, 2.322, 2.500,
                 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191), dtype=float)
    bins["beta"]["n"] = len(bins["beta"]["bins"]) - 1
    
    bins["bpt"] = {}
    #   double vx[] = {15, 20, 25, 30, 35, 40, 50, 60, 75, 90, 110, 130, 175, 230,
                 # 300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750,
                 # 2100, 2500, 3000};
    bins["bpt"]["bins"] = np.array((15, 20, 25, 30, 35, 40, 50, 60, 75, 90, 110, 130, 175, 230,
                 300, 400, 500, 600, 700, 850, 1000, 1200, 1450, 1750,
                 2100, 2500, 3000), dtype=float)
    bins["bpt"]["n"] = len(bins["bpt"]["bins"]) - 1
      # double vptd[] = 
    # {15, 21, 28, 37, 49,
     # 59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575,
     # 638, 737, 846, 967, 1101, 1248,
     # 1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103};
    bins["sami_response_pt"] = {}
    bins["sami_response_pt"]["bins"] = np.array((15, 21, 28, 37, 49,
    59, 86, 110, 132, 170, 204, 236, 279, 302, 373, 460, 575,
    638, 737, 846, 967, 1101, 1248,
    1410, 1588, 1784, 2000, 2238, 2500, 2787, 3103), dtype=float)
    bins["sami_response_pt"]["n"] = len(bins["sami_response_pt"]["bins"]) - 1


    return bins

def file_read_lines(file: str, find_ROOT: bool = False) -> List[str]:
    if find_ROOT:
        with open(file) as f:
            return [line.strip() for line in f.readlines() if line.strip().endswith(".root")]
    else:
        with open(file) as f:
            return [line.strip() for line in f.readlines()]

def read_config_file(config_file: str) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)
    return config

def read_trigger_config(config_file: str) -> Dict:
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(config_file)

    triggers = {}
    for section in config.sections():
        if "filter" not in config[section]:
            triggers[section] = section
        else:
            triggers[section] = "(" + config[section]["filter"] + " && " + section + ")"

    return triggers
