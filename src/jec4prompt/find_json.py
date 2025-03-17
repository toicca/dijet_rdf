import argparse
import json
import os


def update_state(state):
    add_find_json_parser(state.subparsers)
    state.valfuncs["find_json"] = validate_args
    state.commands["find_json"] = run


def add_find_json_parser(subparsers):
    find_json_parser = subparsers.add_parser(
        "find_json", help="Find JSON File appropriate for given run"
    )
    find_json_parser.add_argument(
        "--json_files",
        required=True,
        type=str,
        help="Comma separated \
            list of json files",
    )
    find_json_parser.add_argument(
        "--run_range", required=True, type=str, help="Run number"
    )
    find_json_parser.add_argument("--out", required=False, type=str, help="Output file")


def validate_args(args):
    pass


def run(state):
    args = state.args
    json_files = [s.strip() for s in args.json_files.split(",")]

    run_range = [int(run) for run in args.run_range.split(",")]
    assert len(run_range) == 2

    if args.out:
        output_file = args.out
    else:
        output_file = ""

    state.logger.info(f"json_files: {json_files}")

    newest_run = 0
    newest_json = ""
    partial = True
    for json_file in json_files:
        with open(json_file) as f:
            data = json.load(f)
            runs = [int(r) for r in data.keys()]

            min_run = min(runs)
            max_run = max(runs)

            if (
                run_range[0] >= min_run
                and run_range[1] <= max_run
                and max_run > newest_run
            ):
                newest_run = max_run
                newest_json = json_file

                # JSON file containing the whole run range found.
                # Set partial to False in order to prevent JSON file update
                # on partially matching ranges.
                partial = False
            elif (
                run_range[1] >= min_run
                and run_range[1] <= max_run
                and max_run > newest_run
                and partial
            ):
                # If no JSON file containing the whole run range is found,
                # use JSON file containing the latest run of the given range.
                newest_run = max_run
                newest_json = json_file

    if newest_json != "" and output_file != "":
        os.system(f"cp -r {newest_json} {output_file}")
    elif newest_json != "":
        output_file = newest_json.split("/")[-1]
        os.system(f"cp -r {newest_json} {output_file}")

    state.logger.info(f"Newest JSON file: {output_file}")
    state.logger.info(f"Newest run: {newest_run}")
    state.logger.info(f"Run range: {run_range}")
