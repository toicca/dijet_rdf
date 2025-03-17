#!/usr/bin/env python3
import os
from collections import defaultdict
import argparse


def update_state(state):
    add_find_newest_parser(state.subparsers)
    state.valfuncs["find_newest"] = validate_args
    state.commands["find_newest"] = run


def add_find_newest_parser(subparsers):
    find_newest_parser = subparsers.add_parser(
        "find_newest",
        help="Find newest output file in \
            the subdirectories of given root directory",
    )
    find_newest_parser.add_argument(
        "--root_directory",
        type=str,
        help="Directory to search for \
            files in",
    )
    find_newest_parser.add_argument(
        "--starts_with",
        type=str,
        help="Choose a prefix for the files \
            to search for",
    )
    find_newest_parser.add_argument(
        "--ends_with",
        type=str,
        help="Choose a suffix for the files \
            to search for",
    )
    find_newest_parser.add_argument(
        "--spaces",
        action="store_true",
        help="Use spaces instead of \
            commas to separate the file paths",
    )
    find_newest_parser.add_argument(
        "--max_depth",
        type=int,
        help="Depth of files to search for \
            in the directory tree (default: None)",
    )


def validate_args(args):
    pass


def find_newest_files(root_dir, starts_with, ends_with, max_depth=None):
    newest_files = defaultdict(lambda: {"path": None, "mtime": 0})

    # Walk through the directory tree
    for dirpath, _, filenames in os.walk(root_dir):
        # Calculate the current directory depth relative to root_dir
        current_depth = dirpath[len(root_dir.split(os.path.sep)[0]) :].count(
            os.path.sep
        )

        # Check if we exceed the maximum allowed depth
        if max_depth is not None and current_depth > max_depth:
            # Skip processing files in directories beyond max_depth
            continue

        for filename in filenames:
            if filename.startswith(starts_with):  #  and filename.endswith(ends_with):
                file_path = os.path.join(dirpath, filename)
                mtime = os.path.getmtime(file_path)

                # Check if this file is the newest in its directory
                if mtime > newest_files[dirpath]["mtime"]:
                    newest_files[dirpath]["path"] = file_path
                    newest_files[dirpath]["mtime"] = mtime

    # Return a list of newest file paths
    return [entry["path"] for entry in newest_files.values()]


def run(state):
    args = state.args
    # Specify the root directory to search within
    root_directory = "/eos/user/j/jecpcl/public/test/jec_perIntLumi"
    starts_with = "JEC4PROMPT"
    ends_with = "_plain.root"
    max_depth = None

    if args.root_directory:
        root_directory = args.root_directory
    if args.starts_with:
        starts_with = args.starts_with
    if args.ends_with:
        ends_with = args.ends_with
    if args.max_depth:
        max_depth = args.max_depth
    newest_files = find_newest_files(root_directory, starts_with, ends_with, max_depth)
    if not args.spaces:
        # Print the list comma-separated
        print(",".join(newest_files))
        state.logger.debug(f"Newest files: {newest_files}")
    else:
        # Print the list space-separated
        print(" ".join(newest_files))
        state.logger.debug(f"Newest files: {newest_files}")
