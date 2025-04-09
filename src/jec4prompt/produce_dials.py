from cmsdials import Dials
from cmsdials.filters import (
    RunFilters,
    MLModelsIndexFilters,
)
import json
import argparse

from jec4prompt.utils.processing_utils import file_read_lines, get_rdf, get_files

import logging
logger = logging.getLogger(__name__)

def update_state(state, add_parser=False):
    add_dials_parser(state.subparsers, add_parser=add_parser)
    state.valfuncs["produce_dials"] = validate_args
    state.commands["produce_dials"] = run

def add_dials_parser(subparsers, add_parser=False):
    if add_parser:
        parser = subparsers.add_parser(
            "produce_dials",
            help="Produce DIALS JSON file for a given dataset and run range",
        )
        parser.set_defaults(nThreads=8,
                        progress_bar=True,
                        filepaths='',
                        redirector='root://xrootd-cms.infn.it//',
                        is_local=False)
    else:
        parser = subparsers


    parser.add_argument(
        "--workspace",
        "-ws",
        type=str,
        required=True,
        help="DIALS workspace name",
    )
    parser.add_argument(
        "--json_file",
        type=str,
        default="dials.json",
        help="Output JSON file name",
    )
    parser.add_argument(
        "--default_json",
        type=str,
        default="",
        help="Default JSON file for filtering"
    )
    parser.add_argument(
        "--client_secret",
        type=str,
        default="",
        help="Client secret for DIALS authentication",
    )
    files_parser = parser.add_mutually_exclusive_group(required=True)
    files_parser.add_argument(
        "--filelist",
        type=str,
        help="Comma separated list of input files to parse dataset information",
    )
    files_parser.add_argument(
        "--filepaths",
        type=str,
        help="Path to a .txt file listing the input files"
    )
    parser.add_argument(
        "--nThreads",
        type=int,
        default=8,
        help="Number of threads to use for processing",
    )
    parser.add_argument(
        "--progress_bar",
        action="store_true",
        default=True,
        help="Show progress bar",
    )


def validate_args(args):
    pass

def find_runs(rdf):
    runs = set(rdf.AsNumpy(["run"])["run"])
    return runs

def full_json(rdf):
    runs_ls = rdf.AsNumpy(["run", "luminosityBlock"])

    # For each run find the first and last luminosity block
    run_dict = {}
    for run, ls in zip(runs_ls["run"], runs_ls["luminosityBlock"]):
        if run not in run_dict:
            run_dict[f'{run}'] = [[int(ls), int(ls)]]
        else:
            if ls < run_dict[f'{run}'][0][0]:
                run_dict[f'{run}'][0][0] = int(ls)
            if ls > run_dict[f'{run}'][0][1]:
                run_dict[f'{run}'][0][1] = int(ls)

    return run_dict

def get_datasets(das_files):
    datasets = set()
    regs = set()
    for file in das_files:
        if "/store/" in file:
            base = file.split("/store/data/")[1]
            parts = base.split("/")
            era = parts[0]
            dataset = parts[1]
            datatype = parts[2]
            secondary = parts[3]

            das_dataset = f"/{dataset}/{era}-{secondary}/DQMIO"# {datatype}"
            datasets.add(das_dataset)
            reg = f"/{dataset[:-1]}.*/{era}-{secondary}/DQMIO"
            regs.add(reg)

    return datasets, list(regs)
        

def run(args):
    # Check if args is argparse.Namespace
    if not isinstance(args, argparse.Namespace):
        args = args.args

    # Initialize DIALS
    logger.info("Initializing DIALS credentials")
    if args.client_secret:
        from cmsdials.auth.secret_key import Credentials
        logger.info("Using client secret for DIALS authentication")
        creds = Credentials(token=args.client_secret)
    else:
        from cmsdials.auth.bearer import Credentials
        logger.info("Using interactive DIALS authentication")
        creds = Credentials.from_creds_file()
    dials = Dials(creds, workspace=args.workspace.lower())

    # Find information about the input files
    logger.info("Loading RDF")
    filelist = get_files(args)
    rdf = get_rdf(args)
    logger.info("Finding runs in files")
    runs_in = find_runs(rdf) # Runs in files
    logger.info("Finding datasets in files")
    datasets, _ = get_datasets(filelist) # Datasets in files and regexs

    if not args.default_json:
        logger.info("Finding all lumisections")
        full_json_dict = full_json(rdf) # All possible runs and luminosity blocks
    else:
        logger.info(f"Using {args.default_json} as backup JSON")
        full_json_dict = {}
        with open(args.default_json, "r") as f:
            full_json_dict = json.load(f)
    
    all_runs = set()
    dataset_ids = set()
    for dataset in datasets:
        logger.info(f"Finding DIALS runs in dataset {dataset}")
        # Get the runs in the dataset
        runs = dials.run.list_all(
            RunFilters(page_size=500, dataset=dataset)
        ).results
        logger.debug(f"Found {len(runs)} runs for dataset {dataset}")

        for run in runs:
            all_runs.add(run.run_number)
            dataset_ids.add(run.dataset_id)

    logger.debug(f"All runs in DIALS for datasets {datasets}:\n{all_runs}")

    runs_to_analyze = all_runs.intersection(runs_in, all_runs)
    runs_not_in_dials = runs_in.difference(runs_to_analyze)

    logger.debug(f"Runs to analyze:\n{runs_to_analyze}")
    logger.debug(f"Runs not in DIALS\n{runs_not_in_dials}")

    models_in_dials = dials.ml_models_index.list_all(
        MLModelsIndexFilters(page_size=500, active=True)
    )
    models_in_dials = models_in_dials.to_pandas()
    logger.debug(f"Models in DIALS:\n{models_in_dials}")

    if models_in_dials.empty:
        logger.warning("No models in DIALS, exiting")
        logger.warning(f"Dumping JSON file with all LS as {args.json_file}")

        with open(args.json_file, "w") as f:
            json.dump(full_json_dict, f)
        return

    # Create the JSON file
    mljson = dials.ml_bad_lumis.golden_json(
        model_id__in=models_in_dials.model_id.tolist(),
        dataset_id__in=list(dataset_ids),
        run_number__in=list(runs_to_analyze),
    )

    # If a run is not in the DIALS database, add it to the JSON file
    for run in runs_not_in_dials:
        logger.warning(f"Run {run} not in DIALS")
        if f"{run}" in full_json_dict:
            logger.info(f"Adding run {run} to ML json from default json")
            mljson[f"{run}"] = full_json_dict[f"{run}"]
        else:
            logger.info(f"Run {run} not found in DIALS or default json, skipping it.")

    logger.info(f"Dumping DIALS json to {args.json_file}")
    with open(args.json_file, "w") as f:
        json.dump(mljson, f)


if __name__ == "__main__":
    import sys
    arg_parser = argparse.ArgumentParser(
        description="Produce DIALS JSON file for a given dataset and run range",
    )
    add_dials_parser(arg_parser)
    args = arg_parser.parse_args()
    sys.exit(run(args))




