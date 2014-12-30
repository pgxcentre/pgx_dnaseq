#!/usr/bin/env python
"""The main script that do everything."""


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault", "Valérie Normand"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


import os
import sys
import logging
import argparse
from collections import defaultdict

import fgzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pgx_dna_seq import ProgramError


def main():
    """The main function."""
    # Creating the option parser
    desc = ("Produces a pretty read base quality distribution "
            "plot (version {}).".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    try:
        # Parsing the options
        args = parse_args(parser)
        check_args(args)

        # Adding the logging capability
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=[logging.StreamHandler(),
                      logging.FileHandler(args.log, mode="w")]
        )
        logging.info("Logging everything into '{}'".format(args.log))

        # Reading the two FASTQ files
        quality = read_quality_from_fastq(*args.input_files)

        # Computing the quantiles
        quantiles = compute_percentiles(quality)

        plot_quantiles(quantiles, args.output)

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        logging.warning("Cancelled by user")
        sys.exit(0)

    # Catching the ProgramError
    except ProgramError as e:
        logging.error(e.message)
        parser.error(e.message)


def plot_quantiles(data, output_filename):
    """Plot quantiles."""
    logging.info("Plotting percentiles")

    # The figure and axe
    figure, axe = plt.subplots(1, 1, figsize=(12, 6))

    # Filling
    axe.fill_between(data.pos, data.q5, data.q95, color="#BFBFBF")
    axe.fill_between(data.pos, data.q25, data.q75, color="#7F7F7F")

    # Creating the boxplots from scratch
    axe.plot(data.pos, data.q50, "-", lw=4, color="#CC0000")

    # Adding the labels
    axe.set_title("(A) Pre-Trimming Quality Distribution", weight="bold",
                  fontsize=16)
    axe.set_xlabel("Base Position", weight="bold")
    axe.set_ylabel("PHRED Score", weight="bold")

    # Setting the PHRED score limit
    axe.set_xlim(0, data.pos.max())
    axe.set_ylim(0, data[["q5", "q25", "q50", "q75", "q95"]].max().max())

    # Adding a legend
    lgp = mpatches.Patch(color="#BFBFBF", label="Q5-Q95")
    dgp = mpatches.Patch(color='#7F7F7F', label="Q25-Q75")
    rp = mpatches.Patch(color='#CC0000', label="Median")
    axe.legend(handles=[lgp, dgp, rp], loc="lower left", ncol=3)

    # Adding the grid
    axe.grid(True)

    # Saving the figure
    plt.savefig(output_filename, figure=figure, bbox_inches="tight")
    plt.close(figure)


def compute_percentiles(values):
    """Computes percentiles."""
    logging.info("Computing percentiles")

    # We need to compute the percentiles for each read position
    final_data = []
    pos = []
    for i in range(len(values)):
        # Reconstructing the data
        phreds = np.array([], dtype=int)
        for phred, count in values[i].items():
            phreds = np.r_[phreds, np.zeros(count, dtype=int) + phred]

        # Getting the percentiles
        percentiles = get_percentiles(phreds, [5, 25, 50, 75, 95])

        # Saving the data
        pos.append(i + 1)
        final_data.append(percentiles)

    # Creating the data frame
    final_data = pd.DataFrame(final_data,
                              columns=["q5", "q25", "q50", "q75", "q95"])
    final_data["pos"] = pos

    return final_data


def get_percentiles(data, percentiles):
    """Returns the percentiles of a weighted dataset."""
    return [np.percentile(data, i) for i in percentiles]


def read_quality_from_fastq(filename1, filename2):
    """Reads quality line from two files."""
    open_func = open
    if filename1.endswith(".gz"):
        open_func = fgzip.open

    # The final data
    final_data = [defaultdict(int) for i in range(600)]

    nb_reads = 0
    max_length = -9
    logging.info("Reading {}".format(filename1))
    with open_func(filename1) as i_file:
        # Reading each 4 lines
        for i, line in enumerate(i_file):
            if i % 4 != 3:
                continue

            nb_reads += 1

            qual = [ord(i) - 33 for i in line.decode().rstrip("\n")]
            for j, value in enumerate(qual):
                final_data[j][value] += 1

            if len(qual) > max_length:
                max_length = len(qual)

    open_func = open
    if filename2.endswith(".gz"):
        open_func = fgzip.open

    logging.info("Reading {}".format(filename2))
    with open_func(filename2) as i_file:
        # Reading each 4 lines
        for i, line in enumerate(i_file):
            if i % 4 != 3:
                continue

            nb_reads += 1

            qual = [ord(i) - 33 for i in line.decode().rstrip("\n")]
            for j, value in enumerate(qual):
                final_data[j][value] += 1

            if len(qual) > max_length:
                max_length = len(qual)

    logging.debug("Max length = {}".format(max_length))
    logging.info("Read a total of {:,d} reads".format(nb_reads))

    return final_data[:max_length]


def check_args(args):
    """Checks the arguments and options."""
    # Checking the files
    for filename in args.input_files:
        if not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))
        if (not filename.endswith(".fastq")) and (not filename.endswith(".fastq.gz")):
            raise ProgramError("{}: not a FASTQ file".format(filename))


def parse_args(parser):
    """Parses the command line options and arguments."""
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__))
    parser.add_argument("--debug", action="store_true",
                        help="Set the logging to debug")
    parser.add_argument("--log", type=str, metavar="LOGFILE",
                        default="read_quality_graph.log",
                        help="The log file [%(default)s]")

    # The input files
    group = parser.add_argument_group("Input Files")
    group.add_argument("-i", "--input", type=str, metavar="FILE", nargs=2,
                       dest="input_files",
                       help="The input FASTQ or FASTQ.GZ files")

    # The result file
    group = parser.add_argument_group("Result File")
    group.add_argument("-o", "--output", type=str, metavar="FILE",
                       default="read_quality.pdf",
                       help="The name of the output file [%(default)s]")

    return parser.parse_args()


if __name__ == "__main__":
    main()

