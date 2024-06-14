# ----------------------------------------------------------------------------
# Copyright (c) 2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
from typing import List

from bs4 import BeautifulSoup as BS

EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def get_full_path(file_name):
    """
    Convert a file name or a relative path to its absolute path.

    Args:
    file_name (str): The name of the file or the relative path to the file.

    Returns:
    str: The absolute path of the file.
    """

    return os.path.abspath(file_name)


def create_directory(path):
    """
    Create a directory at the specified path. If the directory already exists,
    no action is taken.

    Args:
    path (str): The path where the directory is to be created.

    Returns:
    bool: True if the directory was created, False if it already existed.
    """

    try:
        os.makedirs(path)
        return True
    except FileExistsError:
        return False


def run_command(cmd, verbose=True):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")
    subprocess.run(cmd, check=True)


def run_commands_with_pipe(cmd1, cmd2, cmd3, outfile_path, verbose=True):
    """Runs two consecutive commands using a pipe"""
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(f'{" ".join(cmd1)} | {" ".join(cmd2)} | {" ".join(cmd3)}', end="\n\n")

    # Start the first command
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    # Pipe the output of cmd1 to cmd2
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits

    # Pipe the output of cmd2 to cmd3, but this time manage the output file ourselves
    with open(outfile_path, "wb") as outfile:
        p3 = subprocess.Popen(cmd3, stdin=p2.stdout, stdout=outfile)
        p2.stdout.close()  # Allow p2 to receive a SIGPIPE if p3 exits

    # Wait for the processes to complete
    p1.wait()
    p2.wait()
    p3.wait()


def _construct_param(arg_name):
    """Converts argument name into a command line parameter."""
    return f'--{arg_name.replace("_", "-")}'


def _process_common_input_params(processing_func, params: dict) -> List[str]:
    """Converts provided arguments and their values.

    Conversion is entirely dependent on the passed 'processing_func'
    that processes individual arguments. The output is a list of
    parameters with their values that can be directly passed to the
    respective command.

    Arguments without any value are skipped.
    Arguments of boolean type are only appended if True.
    Any other argument is processed using the 'processing_func' and
    appended to the final list.

    Args:
        processing_func: Function to be used for formatting a single argument.
        params (dict): Dictionary of parameter: value pairs to be processed.

    Returns:
        processed_args (list): List of processed arguments and their values.

    """
    processed_args = []
    for arg_key, arg_val in params.items():
        if not arg_val:
            continue
        else:
            processed_args.extend(processing_func(arg_key, arg_val))
    return processed_args


def _remove_html_element(fp: str, tag: str, elem_id: str):
    """Removes an HTML tag and its contents.

    Uses BeautifulSoup to open an HTML file, find an element by tag and
    its id and remove it, if found. The original file will be overwritten
    by its modified version.

    Args:
         fp (str): Path to the original HTML file.
         tag (str): Type of the tag to be removed.
         elem_id (str): ID of the element (tag) to be removed.
    """
    with open(fp, "r") as r:
        soup = BS(r.read(), "html.parser")
        element = soup.find(tag, id=elem_id)
        if element:
            element.decompose()
    os.remove(fp)

    with open(fp, "w") as r:
        r.write(str(soup))


def _modify_links(fp: str):
    """Modifies all "a" tags to automatically open in a new tab rather
        than the original iFrame.

    Uses BeautifulSoup to open an HTML file, find all "a" tags and
    add a "target" property. The original file will be overwritten
    by its modified version.

    Args:
         fp (str): Path to the original HTML file.
    """
    with open(fp, "r") as r:
        soup = BS(r.read(), "html.parser")
        links = soup.find_all("a")
        for line in links:
            line["target"] = "_blank"
    os.remove(fp)

    with open(fp, "w") as r:
        r.write(str(soup))


def _get_sample_from_path(fp):
    """Extracts sample name from a contig's file path."""
    return os.path.basename(fp).rsplit("_contigs.fa", maxsplit=1)[0]
