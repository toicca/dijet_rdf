import argparse
import json
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot producer for dijet_rdf: https://github.com/toicca/dijet_rdf")

    parser.add_argument("--json_files", required=True, type=str, help="Comma separated list of json files")

    parser.add_argument("--run", required=True, type=str, help="Run number")

    parser.add_argument("--out", required=True, type=str, help="Output file")

    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = parse_arguments()
    json_files = [s.strip() for s in args.json_files.split(",")]
    run = int(args.run)
    output_file = args.out

    print(f"json_files: {json_files}")
    print(f"run: {run}")
    print(f"out: {output_file}")

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

    if newest_json != "":
        os.system(f"cp -r {newest_json} {output_file}")