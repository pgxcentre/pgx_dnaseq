#!/usr/bin/env python3

"""coverage_graph: A module to plot NGS coverage."""

__author__ = "Louis-Philippe Lemieux Perreault and Abdellatif Daghrach"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault", "Abdellatif Daghrach",
               "Michal Blazejczyk"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"

import os
import re
import sys
import argparse
from csv import QUOTE_NONE
from subprocess import Popen, PIPE

import numpy as np
import pandas as pd


# The version of the script
prog_version = "0.1"

def main():
    """The main function."""
    # Getting and checking the options
    args = parse_args()
    check_args(args)

    # Reading the bed file to get the number of nucleotide that should be
    # covered
    nb_bases = read_bed(args.bed)

    # Getting the depth
    sample_depth = None
    if args.bam is not None:
        sample_depth = compute_sample_depth(args)
    else:
        sample_depth = read_depth(args.depth_file)

    # Plots the graph
    sample_list = list(sample_depth.columns)
    plot_depth(sample_depth, nb_bases, sample_list, args)

def read_bed(filename):
    """Reads a BED file and computes the number of nucleotide included."""
    # Reading the BED file
    bed_file = pd.read_csv(filename, sep="\t", names=["chr", "start", "end"],
                           usecols=range(3)).sort(columns=["chr", "start",
                                                           "end"])

    # Checking that the file is merged
    for i in range(1, len(bed_file)):
        if bed_file.chr.iloc[i] != bed_file.chr.iloc[i-1]:
            continue
        if bed_file.start.iloc[i] + 1 <= bed_file.end.iloc[i-1]:
            m = "{}: regions should be merged".format(filename)
            raise ProgramError(m)

    # Computing the length of the regions (no need to add 1, since the starting
    # position is in 0 format)
    bed_file["length"] = bed_file.end - bed_file.start

    # Returns the number of nucleotides in the targeted regions
    return bed_file.length.sum()

def plot_depth(depth, nb_bases, samples, options):
    # Importing
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    plt.ioff()

    """Plots the depth for all samples."""
    # The figure and axe
    fig, ax = plt.subplots(1, 1, figsize=(14, 8.5))

    # Adding the grids
    ax.grid(color='#8E8E8E', linestyle=':', linewidth=1, which="major")
    ax.set_axisbelow(True)

    # The labels
    ax.set_title(("Sample Depth (MapQ={mapq}, BaseQ={baseq} "
                  "MaxDepth={bam_depth})".format(**vars(options))))
    ax.set_xlabel("Read Depth")
    ax.set_ylabel("Base Proportion")

    # Setting the X limit
    if options.max_depth is not None:
        ax.set_xlim(0, options.max_depth)

    # Plotting for each sample
    for sample in samples:
        # Computing the bins and the cumulative
        bins = np.bincount(depth[sample])
        cumul = np.array([np.sum(bins[i:]) for i in range(len(bins))],
                         dtype=int) / nb_bases

        # Plotting
        plt.plot(np.arange(len(cumul)), cumul, lw=2,
                 label=os.path.basename(sample.split(".")[0]))

    # Plotting the legend
    ax.legend(loc="best",fancybox=True, ncol=3, shadow=True)

    fig.savefig("{}.png".format(options.out), bbox_inches="tight", dpi=300)
    plt.close(fig)

    # Saving the partial data (if required)
    if options.depth_file is None:
        depth.to_csv("{}.txt".format(options.out), sep="\t", index=False,
                    encoding="ascii", quoting=QUOTE_NONE)

def read_depth(filename):
    """Reads the depth from a pre-computed file (from this script)."""
    depth = pd.read_csv(filename, sep="\t", quoting=QUOTE_NONE,
                        encoding="ascii")

    return depth

def compute_sample_depth(options):
    """Computes the sample depth using samtools (and extract relevant info."""
    # The number of sample
    samples = options.bam
    nb_samples = len(samples)

    # The command
    command = ["samtools", "mpileup", "-A", "-d", str(options.bam_depth), "-q",
               str(options.mapq), "-Q", str(options.baseq)]

    # If there is a bed, we add it
    if options.bed is not None:
        command.extend(["-l", options.bed])

    # Adds the input files
    command.extend(options.bam)

    # Launching the subprocess
    p = Popen(command, stdout=PIPE)

    # Reading the output from samtools
    depth = pd.read_csv(p.stdout, sep="\t", names=samples, quoting=QUOTE_NONE,
                        encoding="ascii",
                        usecols=range(3, (nb_samples * 3) + 3, 3))

    # Closing the PIPE
    p.stdout.close()

    return depth

def check_args(args):
    """Checks the arguments and options.

    :param args: an object containing the options and arguments of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    # Checking the BAM files
    if (args.bam is None) and (args.depth_file is None):
        m = "use one of --bam or --depth"
        raise ProgramError(m)

    elif (args.bam is not None) and (args.depth_file is not None):
        m = "use either --bam or --depth"
        raise ProgramError(m)

    if args.bam is not None:
        for filename in args.bam:
            # Is it a BAM or a SAM?
            if re.fullmatch(".*\.[bs]am$", filename) is None:
                m = "{}: not a bam file".format(filename)
                raise ProgramError(m)

            # Does the file exist
            if not os.path.isfile(filename):
                m = "{}: no such file".format(filename)
                raise ProgramError(m)

    elif args.depth_file is not None:
        if not os.path.isfile(args.depth_file):
            m = "{}: no such file".format(args.depth_file)
            raise ProgramError(m)

    # Checking the BED file
    if not args.bed.endswith(".bed"):
        m = "{}: no a bed format".format(args.bed)
        raise ProgramError(m)

    if not os.path.isfile(args.bed):
        m = "{}: no such file".format(args.bed)
        raise ProgramError(m)

    # Checking the qualities
    if args.mapq < 0:
        m = "{}: invalid map quality".format(args.mapq)
        raise ProgramError(m)
    if args.baseq < 0 or args.baseq > 41:
        m = "{}: invalid base quality".format(args.baseq)
        raise ProgramError(m)

    # Checking that the max depth is higher than 0
    if args.max_depth is not None:
        if args.max_depth <= 0:
            m = "{}: invalid maximal depth".format(args.max_depth)
            raise ProgramError(m)
            
    return True

def parse_args():
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the different
              options.

    ================   =======  ================================================
        Options         Type                      Description
    ================   =======  ================================================
    ``--depth-file``   string   Results from this script (to redo the plot
                                faster)
    ``--bam``          string   Input BAM file(s) (one or more, separated by
                                spaces)
    ``--bed``          string   BED file to restrict to targeted regions
    ``-q``             int      skip alignments with mapQ smaller than INT
    ``-Q``             int      skip bases with baseQ/BAQ smaller than INT
    ``-d``             int      max per-BAM depth to avoid excessive memory
                                usage
    ``--max-depth``    int      The maximal depth to plot (in order to zoom in
                                the plots)
    ``--out``          string   The prefix of the output file
    ================   =======  ================================================

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    return parser.parse_args()


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.
    
    :param msg: the message to print to the user before exiting.

    :type msg: string

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user.

        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        return self.message


# The parser object
desc = "Plots NGS coverage (version {}).".format(prog_version)
parser = argparse.ArgumentParser(description=desc)

parser.add_argument("--version", action="version",
                    version="%(prog)s {}".format(prog_version))

# The input files
group = parser.add_argument_group("Input Files")
group.add_argument("--depth-file", type=str, metavar="FILE",
                   help="Results from this script (to redo the plot faster)")
group.add_argument("--bam", type=str, metavar="BAM", nargs="+",
                   help="Input BAM file(s) (one or more, separated by spaces)")
group.add_argument("--bed", type=str, metavar="BED", required=True,
                   help="BED file to restrict to targeted regions")

# The MPILEUP options
group = parser.add_argument_group("MPILEUP Options")
group.add_argument("-q", type=int, metavar="INT", default=0, dest="mapq",
                   help=("skip alignments with mapQ smaller than INT "
                         "[%(default)d]"))
group.add_argument("-Q", type=int, metavar="INT", default=13, dest="baseq",
                   help=("skip bases with baseQ/BAQ smaller than INT "
                         "[%(default)d]"))
group.add_argument("-d", type=int, metavar="INT", default=250, dest="bam_depth",
                   help=("max per-BAM depth to avoid excessive memory usage "
                         "[%(default)d]"))

# Plotting options
group = parser.add_argument_group("Plotting Options")
group.add_argument("--max-depth", type=int, metavar="INT",
                   help=("The maximal depth to plot (in order to zoom in the "
                         "plots) [None]"))

# The output
group = parser.add_argument_group("Output Options")
group.add_argument("-o", "--out", metavar="FILE", default="depth",
                   help="The name of the output file [%(default)s]")

# Calling the main, if necessary
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
