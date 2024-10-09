import ROOT
import configparser
import os
import subprocess
import json
# import tomllib
from processing_utils import find_site, get_bins

def create_histogram(rdf, hist_config, bins):
    if hist_config["type"] == "Histo1D":
        return rdf.Histo1D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"]),
                            hist_config["x_val"], "weight")
    elif hist_config["type"] == "Histo2D":
        return rdf.Histo2D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], "weight")
    elif hist_config["type"] == "Histo3D":
        return rdf.Histo3D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"],
                            bins[hist_config["z_bins"]]["n"], bins[hist_config["z_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], hist_config["z_val"], "weight")
    elif hist_config["type"] == "Profile1D":
        return rdf.Profile1D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], "weight")
    elif hist_config["type"] == "Profile2D":
        return rdf.Profile2D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], hist_config["z_val"], "weight")
    elif hist_config["type"] == "Profile3D":
        return rdf.Profile3D((hist_config["name"], hist_config["title"],
                            bins[hist_config["x_bins"]]["n"], bins[hist_config["x_bins"]]["bins"],
                            bins[hist_config["y_bins"]]["n"], bins[hist_config["y_bins"]]["bins"],
                            bins[hist_config["z_bins"]]["n"], bins[hist_config["z_bins"]]["bins"]),
                            hist_config["x_val"], hist_config["y_val"], hist_config["z_val"], "weight")
    else:
        raise ValueError(f"Unknown histogram type: {hist_config['type']}")

def make_histograms(config):
    ROOT.EnableImplicitMT(config['nThreads'])

    events_chain = ROOT.TChain("Events")
    runs_chain = ROOT.TChain("Runs")

    # Load the files
    for file in config['filelist']:
        if not config['is_local']:
            # Find the available T1, T2 sites
            site_paths = find_site(file)
            if not site_paths:
                print(f"Failed to find site for file '{file}'")
                continue
            for site in site_paths:
                # Test if the file is accessible
                try:
                    f = ROOT.TFile.Open(site_paths[site], "READ")
                    f.Close()
                    events_chain.Add(site_paths[site])
                    runs_chain.Add(site_paths[site])
                    break
                except:
                    continue
        else:
            events_chain.Add(file)
            runs_chain.Add(file)

    events_rdf = ROOT.RDataFrame(events_chain)
    runs_rdf = ROOT.RDataFrame(runs_chain)

    if config['progress_bar']:
        ROOT.RDF.Experimental.AddProgressBar(events_rdf)

    # with open(config['histogram_config'], 'rb') as f:
        # hist_config = tomllib.load(f)
    hist_config = configparser.ConfigParser()
    hist_config.read('histograms.ini')
    hist_config = dict(hist_config)
    
    bins = get_bins()

    histograms = {}

    for hist in hist_config:
        if hist.lower() == "default":
            continue
        histograms[hist] = create_histogram(events_rdf, hist_config[hist], bins).GetValue()

    return histograms

def get_values(histograms):
    values = {}
    for hist in histograms:
        values[hist] = histograms[hist].GetValue()
    return values

def save_histograms(histograms, config):
    output_file = ROOT.TFile(f"J4PHists_runs{config['run_range'][0]}to{config['run_range'][1]}_{config['run_tag']}.root", "RECREATE")

    for hist in histograms:
        histograms[hist].Write()

    output_file.Close()


if __name__ == "__main__":
    config = {
        'filelist': ['J4PSkim_runs379413to379415_20240924.root'],
        'is_local': True,
        'progress_bar': True,
        'histogram_config': 'histograms.toml',
        'run_range': (379413, 379415),
        'run_tag': '20240920',
        'nThreads': 8
    }

    histograms = make_histograms(config)
    save_histograms(histograms, config)


    