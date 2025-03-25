import json
from typing import List

import ROOT
from jec4prompt.utils.processing_utils import file_read_lines, get_bins, read_config_file

response_histos = (
    ("multijet", "MPF", "MPF_multijet_PtRecoilVsEtaVsResponse"),
    ("multijet", "MPF", "MPF_multijet_PtAvpVsEtaVsAvgResponse"),
    ("multijet", "MPF", "MPF_multijet_PtLeadVsEtaVsLeadResponse"),
    ("multijet", "MPF", "MPF_multijet_PtRecoilVsEtaVsRecoilResponse"),
    ("multijet", "DB", "DB_multijet_PtRecoilVsEtaVsResponse"),
    ("multijet", "DB", "DB_multijet_PtRecoilVsEtaVsRecoilResponse"),
    ("multijet", "DB", "DB_multijet_PtAvpVsEtaVsAvgResponse"),
    ("multijet", "DB", "DB_multijet_PtLeadVsEtaVsLeadResponse"),
    ("dijet", "MPF", "MPF_dijet_PtTagVsEtaVsResponse"),
    ("dijet", "MPF", "MPF_dijet_PtTagVsEtaVsTagResponse"),
    ("dijet", "MPF", "MPF_dijet_PtAvgVsEtaVsAvgResponse"),
    ("dijet", "MPF", "MPF_dijet_PtProbeVsEtaVsProbeResponse"),
    ("dijet", "DB", "DB_dijet_PtTagVsEtaVsResponse"),
    ("dijet", "DB", "DB_dijet_PtTagVsEtaVsTagResponse"),
    ("dijet", "DB", "DB_dijet_PtAvgVsEtaVsAvgResponse"),
    ("dijet", "DB", "DB_dijet_PtProbeVsEtaVsProbeResponse"),
)
resolution_histos = (
    ("dijet", "DB", "DB_dijet_PtTagVsEtaVsA"),
    ("dijet", "DB", "DB_dijet_PtProbeVsEtaVsA"),
    ("dijet", "DB", "DB_dijet_PtAvgVsEtaVsA"),
)

derived_histos = ()
# derived_histos = (("multijet", "MPF", "MPF_multijet_PtAvgVsEtaVsB"),
# ("multijet", "MPF", "MPF_multijet_PtRecoilVsEtaVsB"),
# ("multijet", "MPF", "MPF_multijet_PtLeadVsEtaVsB"),
# ("multijet", "DB", "DB_multijet_PtAvgVsEtaVsA"),
# ("multijet", "DB", "DB_multijet_PtRecoilVsEtaVsA"),
# ("multijet", "DB", "DB_multijet_PtLeadVsEtaVsA"),
# )


def update_state(state):
    add_produce_responses_parser(state.subparsers)
    state.valfuncs["produce_responses"] = validate_args
    state.commands["produce_responses"] = run


def add_produce_responses_parser(subparsers):
    responses_parser = subparsers.add_parser(
        "produce_responses",
        help="Produce responses \
            for files produced by JEC4PROMPT analysis",
    )
    responses_files = responses_parser.add_mutually_exclusive_group(required=True)
    responses_files.add_argument(
        "--filelist",
        type=str,
        help="Comma separated list of root files \
            produced by dijet_rdf",
    )
    responses_files.add_argument(
        "-fp",
        "--filepaths",
        type=str,
        help="Comma separated list of \
            text files containing input files (one input file per line).",
    )
    responses_parser.add_argument(
        "-tf",
        "--triggerfile",
        type=str,
        help="Path to the .json file containing \
            triggers.",
    )
    responses_parser.add_argument("--out", type=str, default="", help="Output path")
    responses_parser.add_argument(
        "--config", type=str, default="", help="Path to config file"
    )


def validate_args(args):
    pass


def produce_resolutions(file: str, trigger_list: List[str], output_path: str):
    """
    Resolution producer for dijet_rdf.
    """

    bins = get_bins()

    file = ROOT.TFile(file, "UPDATE")
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers in the file.")
        trigger_keys = file.GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]

    for trg in trigger_list:
        for system, method, histogram in resolution_histos:
            path = f"{trg}/{system}/{method}/"
            resolution_path = f"{trg}/{system}/Resolutions"

            if not file.GetDirectory(resolution_path):
                file.mkdir(resolution_path)

            h = file.Get(path + histogram)
            if not h:
                print(f"Could not find {path + histogram}")
                continue
            # h.RebinY(4)
            xlabel = ""
            if "PtTag" in histogram:
                xlabel = "p_{T, tag}"
            elif "PtProbe" in histogram:
                xlabel = "p_{T, probe}"
            elif "PtAvg" in histogram:
                xlabel = "p_{T, avg}"

            resolutions = ROOT.TH2D(
                "resolutions_"
                + histogram.replace("VsA", "")
                .replace("VsEta", "")
                .replace(f"_{method}_{system}", ""),
                h.GetTitle() + ";" + xlabel + ";#eta_{probe}",
                bins["pt"]["n"],
                bins["pt"]["bins"],
                h.GetNbinsY(),
                h.GetYaxis().GetXmin(),
                h.GetYaxis().GetXmax(),
            )
            projections = ROOT.TH2D(
                "projections_"
                + histogram.replace("VsA", "")
                .replace("VsEta", "")
                .replace(f"_{method}_{system}", ""),
                h.GetTitle(),
                bins["pt"]["n"],
                bins["pt"]["bins"],
                bins["pt"]["n"],
                bins["pt"]["bins"],
            )

            # TH3.FitSlicesZ?
            for i in range(1, h.GetNbinsX() + 1):
                for j in range(1, h.GetNbinsY() + 1):
                    # Fit to the asymmetry in the bin
                    proj = h.ProjectionZ(f"proj_{i}_{j}", i, i, j, j)
                    # Find the mean
                    mean = proj.GetMean()
                    # Find the sigma
                    sigma = 0
                    for k in range(1, proj.GetNbinsX() + 1):
                        sigma += (proj.GetBinContent(k) - mean) ** 2
                    sigma = proj.GetRMS()

                    if sigma > 0:
                        gaussian = ROOT.TF1(
                            "my_gaus", "gaus", mean - sigma, mean + sigma
                        )
                        proj.Fit(gaussian, "Q L N")
                        fsigma = gaussian.GetParameter(2)
                        resolutions.SetBinContent(i, j, fsigma)
                        resolutions.SetBinError(i, j, gaussian.GetParError(2))
                    else:
                        resolutions.SetBinContent(i, j, 0)
                        resolutions.SetBinError(i, j, 0)
                    # Save the projection
                    for k in range(1, proj.GetNbinsX() + 1):
                        projections.SetBinContent(i, k, proj.GetBinContent(k))

            resolutions.SetName(
                method
                + "_projected_resolution_"
                + histogram.replace("VsA", "")
                .replace("VsEta", "")
                .replace(f"_{method}_{system}", "")
            )

            # Save
            file.cd(resolution_path)
            # h2.Write()
            resolutions.Write()
            file.cd()


def produce_responses(file: str, trigger_list: List[str], output_path: str):
    """
    Response producer for dijet_rdf.
    """

    bins = get_bins()

    file = ROOT.TFile(file, "UPDATE")
    if len(trigger_list) == 0:
        print("No triggers provided. Using all triggers in the file.")
        trigger_keys = file.GetListOfKeys()
        trigger_list = [tkey.GetName() for tkey in trigger_keys]

    for trg in trigger_list:
        for system, method, histogram in response_histos:
            path = f"{trg}/{system}/{method}/"
            response_path = f"{trg}/{system}/Responses"

            if not file.GetDirectory(response_path):
                file.mkdir(response_path)

            h = file.Get(path + histogram)
            if not h:
                print(f"Could not find {path + histogram}")
                continue
            # Get the projection w.r.t. y-axis
            h.GetYaxis().SetRangeUser(-2.5, 2.5)
            h2 = h.Project3D("zx")

            # Profile
            h3 = h2.ProfileX()
            h3.SetName(
                method
                + "_projected_response_"
                + histogram.replace("VsResponse", "")
                .replace("VsEta", "")
                .replace(f"_{method}_{system}", "")
            )

            # Save
            file.cd(response_path)
            # h2.Write()
            h3.Write()
            file.cd()

        for system, method, histogram in derived_histos:
            path = f"{trg}/{system}/{method}/"
            response_path = f"{trg}/{system}/Responses"
            if not file.GetDirectory(path):
                file.mkdir(path)
            if not file.GetDirectory(response_path):
                file.mkdir(response_path)

            h = file.Get(path + histogram)
            # Get the projection w.r.t. y-axis
            h.GetYaxis().SetRangeUser(-2.5, 2.5)
            h2 = h.Project3D("zx")

            # # Profile
            h3 = h2.ProfileX().ProjectionX()

            # # Unit histogram
            unit = h3.Clone()
            unit.Reset()
            for i in range(1, unit.GetNbinsX() + 1):
                unit.SetBinContent(i, 1)
                unit.SetBinError(i, 0)

            nominator = unit.Clone()
            nominator.Add(unit, h3, c1=1.0, c2=1.0)
            denominator = unit.Clone()
            denominator.Add(unit, h3, c1=1.0, c2=-1.0)

            h3.Divide(nominator, denominator, 1, 1)

            h3.SetName(
                method
                + "_derived_response_"
                + histogram.replace("VsResponse", "")
                .replace("VsEta", "")
                .replace(f"_{method}_{system}", "")
            )

            # Save
            file.cd(response_path)
            h3.Write()


def run(state):
    args = state.args
    files: List[str] = []

    with open(args.triggerfile, "r") as f:
        trigger_list = json.load(f)[args.channel]

    for trigger in trigger_list:
        trigger_list[trigger] = trigger_list[trigger]["cut"]

    # Split the file list and trigger list if they are given as a string
    if args.filelist:
        files = args.filelist.split(",")
    elif args.filepaths:
        paths = [p.strip() for p in args.filepaths.split(",")]
        for path in paths:
            files.extend(file_read_lines(path, find_ROOT=True))
    else:
        raise ValueError("No file list provided")

    output_path = args.out

    if args.config:
        config_file = args.config
        config = read_config_file(config_file)

    for file in files:
        produce_responses(file, trigger_list, output_path)
        produce_resolutions(file, trigger_list, output_path)
