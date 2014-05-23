#!/usr/bin/env python3

import os
import sys
import argparse

import pgx_dna_seq
from pgx_dna_seq import ProgramError


# The version of the script
prog_version = pgx_dna_seq.__version__

def main():
    """The main function."""
    # Getting and checking the options
    args = parse_args()
    check_args(args)

    # Summarize the options used
    print("\n# Command used:")
    print("{} \\".format(sys.argv[0]))
    for option, value in vars(args).items():
        option = option.replace("_", "-")
        print("    --{} {} \\".format(option, value))
    print()

def check_args(args):
    """Checks the arguments and options.

    :param args: an object containing the options and arguments of the program.

    :type args: :py:class:`argparse.Namespace`

    :returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    # Checking that the configuration file exists
    if not os.path.isfile(args.pipeline_config):
        m = "{}: no such file".format(args.pipeline_config)
        raise ProgramError(m)
    if not os.path.isfile(args.tool_config):
        m = "{}: no such file".format(args.tool_config)
        raise ProgramError(m)

    return True

def parse_args():
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the different
              options.

    ===============   =======  ================================================
        Options        Type                      Description
    ===============   =======  ================================================
    ===============   =======  ================================================

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    return parser.parse_args()

# The parser object
desc = "Execute a NGS pipeline (version {}).".format(prog_version)
parser = argparse.ArgumentParser(description=desc)

# The input file
parser.add_argument("--version", action="version",
                    version="%(prog)s {}".format(prog_version))

group = parser.add_argument_group("Pipeline Configuration")
group.add_argument("--pipeline-config", type=str, metavar="FILE",
                   default="pipeline.conf",
                   help=("The pipeline configuration file (default: "
                         "%(default)s)"))
group.add_argument("--tool-config", type=str, metavar="FILE",
                   default="tools.conf",
                   help="The tools configuration file (default: %(default)s)")
group.add_argument("--use-drmaa", action="store_true", default=False,
                   help=("Use DRMAA to launch the tasks instead of running "
                         "them locally (default: %(default)s)"))
group.add_argument("--nb-process", type=int, metavar="INT", default=1,
                   help=("The number of processes for job execution (allow "
                         "enough if '--use-drmaa' option is used since you "
                         "want at least one job per sample) (default: "
                         "%(default)s)"))

# Calling the main, if necessary
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print >>sys.stderr, "Cancelled by user"
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
