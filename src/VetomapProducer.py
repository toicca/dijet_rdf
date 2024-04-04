import ROOT
from RDFHelpers import *

if __name__ == '__main__':
    input_file = "out/dijet_test_22Mar.root"
    trigger_file = "data/testtrigger.txt"
    
    trigger_list = readTriggerList(trigger_file)

    f = ROOT.TFile(input_file, "UPDATE")
    hname = "PFComposition_EtaVsPhiVsProfileNHF_selected"

    for trg in trigger_list:
        if trg == "HLT_ZeroBias":
            continue
        # out_file = f"out/vetomap_tests/{trg}_vetomap.root"
        
        obj_path = f"{trg}/PFComposition/{hname}"
        obj = f.Get(obj_path)
        assert obj, f"Object not found at {obj_path}"
        assert obj.InheritsFrom("TH2D"), f"Object at {obj_path} is not a TH2D or derived class"

        isProf2D = obj.InheritsFrom("TProfile2D")
        h2 = obj.ProjectionXY() if isProf2D else obj
        h2nom = h2.Clone(f"h2nom_{hname}")
        h2abs = h2.Clone(f"h2abs_{hname}")

        # Calculate average JES shift
        if hname == "p2asymm":
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
                        # The doPull variable or its equivalent isn't defined in the snippet provided.
                        # If you have specific logic for doPull, it needs to be adapted here.
                        delta = abs(h2.GetBinContent(i, j) - median) / q68
                        h2abs.SetBinContent(i, j, delta - 1)  # Assuming doPull equivalent is intended
                        h2abs.SetBinError(i, j, h2.GetBinError(i, j) / q68)
                        h2nom.SetBinContent(i, j, delta)
                        h2nom.SetBinError(i, j, h2.GetBinError(i, j) / q68)
            
            htmp.Delete()
            
        # Write the histograms to the output file
        # out_f = ROOT.TFile(out_file, "RECREATE")
        # out_f.cd()
        # Write to the trigger directory
        # TODO
        f.mkdir(trg+"/VetoMap")
        f.cd(trg+"/VetoMap")
        h2nom.Write()
        h2abs.Write()
        f.cd()
        
        

    