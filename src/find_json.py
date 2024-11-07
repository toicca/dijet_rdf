import argparse
import json
import os

def run(args):
    json_files = [s.strip() for s in args.json_files.split(",")]

    run_range = [int(run) for run in args.run_range.split(",")]
    assert(len(run_range) == 2)

    if args.out:
        output_file = args.out
    else:
        output_file = ""

    #print(f"json_files: {json_files}")
    #print(f"run: {run}")

    newest_run = 0
    newest_json = ""
    partial = True
    for json_file in json_files:
        with open(json_file) as f:
            data = json.load(f)
            runs = [int(r) for r in data.keys()]

            min_run = min(runs)
            max_run = max(runs)

            if run_range[0] >= min_run and run_range[1] <= max_run and max_run > newest_run:
                newest_run = max_run
                newest_json = json_file

                # JSON file containing the whole run range found.
                # Set partial to False in order to prevent JSON file update
                # on partially matching ranges.
                partial = False
            elif run_range[1] >= min_run and run_range[1] <= max_run \
                    and max_run > newest_run and partial:
                # If no JSON file containing the whole run range is found,
                # use JSON file containing the latest run of the given range.
                newest_run = max_run
                newest_json = json_file

    if newest_json != "" and output_file != "":
        os.system(f"cp -r {newest_json} {output_file}")
    elif newest_json != "":
        output_file = newest_json.split("/")[-1]
        os.system(f"cp -r {newest_json} {output_file}")

    print(output_file)
