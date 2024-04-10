import ROOT
from RDFHelpers import file_read_lines, read_config_file, get_bins
from typing import List
import argparse, configparser
import numpy as np

hist_names = ("PFComposition_EtaVsPhiVsProfileNEF_selected",
            "PFComposition_EtaVsPhiVsProfileCEF_selected",
            "PFComposition_EtaVsPhiVsProfileCHF_selected",
            "DB_dijet_EtaprobeVsPhiprobeVsAsymmetry",
            "Inclusive_EtaVsPhi_selected"
            )


def parse_arguments():
    parser = argparse.ArgumentParser(description="VetoMap for dijet_rdf: https://github.com/toicca/dijet_rdf")

    files = parser.add_mutually_exclusive_group(required=True)
    files.add_argument("--filelist", type=str, help="Comma separated list of root files produced by dijet_rdf")
    files.add_argument("--filepath", type=str, help="Path to a root file containing a list of output files produced by dijet_rdf")

    triggers = parser.add_mutually_exclusive_group()
    triggers.add_argument("--triggerlist", type=str, help="Comma separated list of triggers for which plots will be produced (default value 'all')")
    triggers.add_argument("--triggerpath", type=str, help="Path to a file containing a list of triggers for which plots will be produced")

    parser.add_argument("--out", type=str, default="", help="Output path")

    parser.add_argument("--config", type=str, default="", help="Path to config file")

    args = parser.parse_args()
    
    return args

def produce_vetomap(input_file: str, trigger_list: List[str], hist_names: List[str], output_path: str):
    """
    VetoMap producer for dijet_rdf.
    Translated from https://github.com/miquork/jecsys3/blob/21fcfb6c4a6fab963650b08019f57d679644ec83/minitools/doJetVetoV2.C
    """
    
    bins = get_bins()
    
    file = ROOT.TFile(input_file, "UPDATE")
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers in the file.")
        trigger_keys = file.GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]
        
    # Get the folder name from hist_names
    methods = [h.split("_")[0] for h in hist_names]

    for trg in trigger_list:
        h2nomsum = ROOT.TH2D("h2nomsum", "Sum of methods", bins["eta"]["n"], bins["eta"]["bins"], bins["phi"]["n"], bins["phi"]["bins"])
        h2abssum = ROOT.TH2D("h2abssum", "Sum of methods", bins["eta"]["n"], bins["eta"]["bins"], bins["phi"]["n"], bins["phi"]["bins"])
        for method, hname in zip(methods, hist_names):            
            obj_path = f"{trg}/{method}/{hname}"
            obj = file.Get(obj_path)
            assert obj, f"Object not found at {obj_path}"
            assert obj.InheritsFrom("TH2D"), f"Object at {obj_path} is not a TH2D or derived class"

            isProf2D = obj.InheritsFrom("TProfile2D")
            h2 = obj.ProjectionXY() if isProf2D else obj
            h2nom = h2.Clone(f"h2nom_{hname}")
            h2abs = h2.Clone(f"h2abs_{hname}")

            # Calculate average JES shift
            if hname == "DB_dijet_EtaprobeVsPhiprobeVsAsymmetry":
                if 'h2jes' not in locals():
                    h2jes = h2.Clone("h2jes")
                else:
                    for i in range(1, h2.GetNbinsX()+1):
                        for j in range(1, h2.GetNbinsY()+1):
                            if h2.GetBinError(i, j) != 0:
                                if h2jes.GetBinError(i, j) != 0:
                                    val1, err1 = h2jes.GetBinContent(i, j), h2jes.GetBinError(i, j)
                                    n1 = 1./pow(err1, 2)
                                    val2, err2 = h2.GetBinContent(i, j), h2.GetBinError(i, j)
                                    n2 = 1./pow(err2, 2)
                                    val = (n1*val1 + n2*val2) / (n1+n2)
                                    err = (n1*err1 + n2*err2) / (n1+n2)
                                    h2jes.SetBinContent(i, j, val)
                                    h2jes.SetBinError(i, j, err)
                                elif h2.GetBinContent(i, j) != 0:
                                    h2jes.SetBinContent(i, j, h2.GetBinContent(i, j))
                                    h2jes.SetBinError(i, j, h2.GetBinError(i, j))

            # Normalize eta strips vs phi
            for i in range(1, h2.GetNbinsX()+1):
                htmp = ROOT.TH1D("htmp", "Distribution of values", 100, -1, -1)
                # htmp.SetBuffer(72)
                n = 0
                for j in range(1, h2.GetNbinsY()+1):
                    if h2.GetBinError(i, j) != 0:
                        htmp.Fill(h2.GetBinContent(i, j))
                        n += 1
                
                if n < 70:
                    for j in range(1, h2.GetNbinsY()+1):
                        h2.SetBinContent(i, j, 0.)
                        h2.SetBinError(i, j, 0.)
                else:
                    probSum = np.array([0.16, 0.50, 0.84], dtype=np.float64)
                    q = np.array([0., 0., 0.], dtype=np.float64)
                    htmp.GetQuantiles(len(probSum), q, probSum)
                    median = q[1]
                    q68 = 0.5*(q[2]-q[0])

                    for j in range(1, h2.GetNbinsY()+1):
                        if h2.GetBinError(i, j) != 0 and q68 != 0:
                            delta = abs(h2.GetBinContent(i, j) - median) / q68
                            h2abs.SetBinContent(i, j, delta - 1)
                            h2abs.SetBinError(i, j, h2.GetBinError(i, j) / q68)
                            h2nom.SetBinContent(i, j, delta)
                            h2nom.SetBinError(i, j, h2.GetBinError(i, j) / q68)
                
                htmp.Delete()
                
            if not file.GetDirectory(trg):
                file.mkdir(trg)
            if not file.GetDirectory(trg+"/VetoMap"):
                file.mkdir(trg+"/VetoMap")
            file.cd(trg+"/VetoMap")
            h2nom.Write()
            h2abs.Write()
            file.cd()
            
    file.Close()
    return

def produce_vetomap2(input_file: str, trigger_list: List[str], hist_names: List[str], output_path: str):
    """
    VetoMap producer for dijet_rdf.
    Translated from https://github.com/miquork/jecsys3/blob/21fcfb6c4a6fab963650b08019f57d679644ec83/minitools/doJetVetoV2.C
    """
    
    bins = get_bins()
    
    file = ROOT.TFile(input_file, "UPDATE")
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers in the file.")
        trigger_keys = file.GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]
        
    # Get the folder name from hist_names
    methods = [h.split("_")[0] for h in hist_names]

    for trg in trigger_list:
        # Create an empty histogram to store the sum of the method histograms
        sum_of_methods = ROOT.TH2D("sum_of_methods", "Sum of methods", bins["eta"]["n"], bins["eta"]["bins"], bins["phi"]["n"], bins["phi"]["bins"])
        veto_map_loose = ROOT.TH2D("veto_map_loose", "Veto map loose", bins["eta"]["n"], bins["eta"]["bins"], bins["phi"]["n"], bins["phi"]["bins"])
        veto_map_medium = ROOT.TH2D("veto_map_medium", "Veto map medium", bins["eta"]["n"], bins["eta"]["bins"], bins["phi"]["n"], bins["phi"]["bins"])
        veto_map_tight = ROOT.TH2D("veto_map_tight", "Veto map tight", bins["eta"]["n"], bins["eta"]["bins"], bins["phi"]["n"], bins["phi"]["bins"])
        spread_of_deltas = ROOT.TH1D("spread_of_deltas", "Spread of deltas", 100, -3, 3)
        deltas = ROOT.TH2D("deltas", "Deltas", bins["eta"]["n"], bins["eta"]["bins"], bins["phi"]["n"], bins["phi"]["bins"])
        
        # Initialize as zeros
        for i in range(1, sum_of_methods.GetNbinsX()+1):
            for j in range(1, sum_of_methods.GetNbinsY()+1):
                sum_of_methods.SetBinContent(i, j, 0)
                sum_of_methods.SetBinError(i, j, 0)
        
        # Get the method histograms
        for method, hname in zip(methods, hist_names):            
            obj_path = f"{trg}/{method}/{hname}"
            obj = file.Get(obj_path)
            assert obj, f"Object not found at {obj_path}"
            assert obj.InheritsFrom("TH2D"), f"Object at {obj_path} is not a TH2D or derived class"

            isProf2D = obj.InheritsFrom("TProfile2D")
            h2 = obj.ProjectionXY() if isProf2D else obj
            h2nom = h2.Clone(f"h2nom_{hname}")
            h2abs = h2.Clone(f"h2abs_{hname}")
            
            
            # Calculate average JES shift
            if hname == "DB_dijet_EtaprobeVsPhiprobeVsAsymmetry":
                for i in range(1, h2.GetNbinsX()+1):
                    for j in range(1, h2.GetNbinsY()+1):
                        if h2.GetBinError(i, j) != 0:
                            if h2.GetBinError(i, j) != 0:
                                val1, err1 = h2.GetBinContent(i, j), h2.GetBinError(i, j)
                                n1 = 1./pow(err1, 2)
                                val2, err2 = h2.GetBinContent(i, j), h2.GetBinError(i, j)
                                n2 = 1./pow(err2, 2)
                                val = (n1*val1 + n2*val2) / (n1+n2)
                                err = (n1*err1 + n2*err2) / (n1+n2)
                                h2.SetBinContent(i, j, val)
                                h2.SetBinError(i, j, err)
                                  
            # Normalize eta strips
            for i in range(1, h2.GetNbinsX()+1):
                # Normalize eta strip to area
                projection = h2.ProjectionY(f"projection_{i}", i, i)
                area = projection.Integral()
                if area != 0:
                    for j in range(1, h2.GetNbinsY()+1):
                        h2.SetBinContent(i, j, h2.GetBinContent(i, j) / area)
                        h2.SetBinError(i, j, h2.GetBinError(i, j) / area)
                
                # Add to the sum of methods
                for j in range(1, h2.GetNbinsY()+1):
                    sum_of_methods.SetBinContent(i, j, sum_of_methods.GetBinContent(i, j) + h2.GetBinContent(i, j))
                    sum_of_methods.SetBinError(i, j, sum_of_methods.GetBinError(i, j) + h2.GetBinError(i, j))
            
        # Normalize the sum of methods to strip mean = 1
        sum_of_methods.Scale(1./len(methods))
                
        # Calculate the veto maps
        
        # If the value of the sum of methods is between mean +- probSum[i]*sigma, set the veto map value to 0, else set it to 1
        # Loop through the eta strips:
        for i in range(1, sum_of_methods.GetNbinsX()+1):
            projection = sum_of_methods.ProjectionY(f"projection_{i}", i, i)
            # Get median by summing bin content and dividing by number of bins
            if projection.GetNbinsX() % 2 == 0:
                median = (projection.GetBinContent(projection.GetNbinsX()//2) + projection.GetBinContent(projection.GetNbinsX()//2 + 1)) / 2
            else:
                median = projection.GetBinContent(projection.GetNbinsX()//2)

            # Get sigma by summing bin errors in quadrature and dividing by number of bins
            sum = 0
            for j in range(1, projection.GetNbinsX()+1):
                sum += pow(projection.GetBinContent(j) - median, 2)
            sigma = np.sqrt(sum / projection.GetNbinsX())
            
            # Check how many sigmas is the point away from the mean
            if sigma != 0:
                for j in range(1, sum_of_methods.GetNbinsY()+1):
                    # Require enough events in the bin
                    if (sum_of_methods.GetBinError(i, j) != 0) and (sum_of_methods.GetBinError(i, j)/median < 0.5):
                        delta = (sum_of_methods.GetBinContent(i, j) - median) / sigma
                        veto_map_loose.SetBinContent(i, j, 0 if ((delta < 5) and (delta > -5)) else 1)
                        veto_map_medium.SetBinContent(i, j, 0 if ((delta < 4) and (delta > -4)) else 1)
                        veto_map_tight.SetBinContent(i, j, 0 if ((delta < 3) and (delta > -3)) else 1)
                        spread_of_deltas.Fill(delta)
                        deltas.SetBinContent(i, j, delta)
                    else:
                        veto_map_loose.SetBinContent(i, j, 0)
                        veto_map_medium.SetBinContent(i, j, 0)
                        veto_map_tight.SetBinContent(i, j, 0)
                        spread_of_deltas.Fill(0)
                        deltas.SetBinContent(i, j, 0)
                    
                    
        
        # Write the veto maps to the file
        if not file.GetDirectory(trg):
            file.mkdir(trg)
        if not file.GetDirectory(trg+"/VetoMap2"):
            file.mkdir(trg+"/VetoMap2")
        file.cd(trg+"/VetoMap2")
        sum_of_methods.Write()
        veto_map_loose.Write()
        veto_map_medium.Write()
        veto_map_tight.Write()
        spread_of_deltas.Write()
        deltas.Write()
        file.cd()
        

        
        


if __name__ == '__main__':
    
    args = parse_arguments()
    
    trigger_list: List[str] = []
    files: List[str] = []
    
    if args.triggerpath:
        trigger_list = file_read_lines(args.triggerpath)
    elif args.triggerlist:
        trigger_list = [s.strip() for s in args.triggerlist.split(',')]

    if args.filepath:
        files = file_read_lines(args.filepath)
    else:
        files = [s.strip() for s in args.filelist.split(',')]

    output_path = args.out

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    for file in files:
        produce_vetomap(file, trigger_list, hist_names, output_path)
        produce_vetomap2(file, trigger_list, hist_names, output_path)
    
        
        

    