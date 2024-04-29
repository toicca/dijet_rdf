#!/usr/bin/env python3
import os
from collections import defaultdict
import argparse

# Specify the root directory to search within
root_directory = '/eos/user/j/jecpcl/public/test/2024C/jec_perIntLumi'
starts_with = 'combined'

def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to find the newest files in subdirectories for dijet_rdf: https://github.com/toicca/dijet_rdf")

    parser.add_argument("--root_directory", required=False, type=str, help="Directory to search for files in")
    parser.add_argument("--starts_with", required=False, type=str, help="Choose a prefix for the files to search for")

    args = parser.parse_args()
    
    return args

def find_newest_files(root_dir, starts_with):
    newest_files = defaultdict(lambda: {"path": None, "mtime": 0})

    # Walk through the directory tree
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith('.root') and filename.startswith(starts_with):
                file_path = os.path.join(dirpath, filename)
                mtime = os.path.getmtime(file_path)

                # Check if this file is the newest in its directory
                if mtime > newest_files[dirpath]["mtime"]:
                    newest_files[dirpath]["path"] = file_path
                    newest_files[dirpath]["mtime"] = mtime

    # Return a list of newest file paths
    return [entry["path"] for entry in newest_files.values()]

if __name__ == "__main__":
    args = parse_arguments()
    if args.root_directory:
        root_directory = args.root_directory
    if args.starts_with:
        starts_with = args.starts_with
    newest_files = find_newest_files(root_directory, starts_with)

    for file_path in newest_files:
        print(file_path)
