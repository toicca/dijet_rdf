import argparse
import json
import os

def run(args):
    json_files = [s.strip() for s in args.json_files.split(",")]
    run = int(args.run)
    if args.out:
        output_file = args.out
    else:
        output_file = ""

    print(f"json_files: {json_files}")
    print(f"run: {run}")

    newest_run = 0
    newest_json = ""
    for json_file in json_files:
        with open(json_file) as f:
            data = json.load(f)
            runs = [int(r) for r in data.keys()]

            min_run = min(runs)
            max_run = max(runs)

            if run >= min_run and run <= max_run and max_run > newest_run:
                newest_run = max_run
                newest_json = json_file

    if newest_json != "" and output_file != "":
        os.system(f"cp -r {newest_json} {output_file}")
    elif newest_json != "":
        output_file = newest_json.split("/")[-1]
        os.system(f"cp -r {newest_json} {output_file}")

    print(f"out: \n{newest_json}")
    # print(f"out: \n{output_file}")
