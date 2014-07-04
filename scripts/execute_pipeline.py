#!/usr/bin/env python3

import os
import re
import sys
import argparse
import __main__
import traceback
from copy import copy

import pgx_dna_seq
from ruffus import originate, formatter, collate, transform, regex
from ruffus import pipeline_printout_graph, pipeline_run
from pgx_dna_seq import ProgramError
from pgx_dna_seq.tool import GenericTool as Tool
from pgx_dna_seq.read_config import read_config_file, get_pipeline_steps


# The version of the script
prog_version = pgx_dna_seq.__version__


def check_input_files(filename):
    """Checks the input file names."""
    split_re = re.compile(r"\s+")
    input_filenames = None
    with open(filename, "r") as i_file:
        input_filenames = [re.split(r"\s+", i.rstrip("\r\n"))
                                                    for i in i_file.readlines()]

    # Checking that all those files exists
    for sample_files in input_filenames:
        for input_filename in sample_files:
            if not os.path.isfile(input_filename):
                m = "{}: no such file".format(input_filename)
                raise ProgramError(m)

    return input_filenames


def rename_func(new_name):
    """Decorator function that renames a function."""
    def decorator(func):
        func.__name__ = new_name
        return func
    return decorator


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

    # Checking that the file containing all input files exists
    if not os.path.isfile(args.input):
        m = "{}: no such file".format(args.input)
        raise ProgramError(m)

    # Checking that the number of process is higher than 0
    if args.nb_process < 1:
        m = "{}: invalid number of process".format(args.nb_process)
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
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s version {}".format(prog_version))

group = parser.add_argument_group("Pipeline Configuration")
group.add_argument("-i", "--input", type=str, metavar="FILE",
                   default="input_files.txt",
                   help=("A file containing the pipeline input files (one "
                         "sample per line, one or more file per sample "
                         "[%(default)s]"))
group.add_argument("-p", "--pipeline-config", type=str, metavar="FILE",
                   default="pipeline.conf",
                   help="The pipeline configuration file [%(default)s]")
group.add_argument("-t", "--tool-config", type=str, metavar="FILE",
                   default="tools.conf",
                   help="The tools configuration file [%(default)s]")
group.add_argument("-d", "--use-drmaa", action="store_true", default=False,
                   help=("Use DRMAA to launch the tasks instead of running "
                         "them locally [%(default)s]"))
group.add_argument("-n", "--nb-process", type=int, metavar="INT", default=1,
                   help=("The number of processes for job execution (allow "
                         "enough if '--use-drmaa' option is used since you "
                         "want at least one job per sample) [%(default)s]"))


# Calling the main, if necessary
if __name__ == "__main__":
    try:
        # Getting and checking the options
        args = parse_args()
        check_args(args)

        # Summarize the options used
        print()
        print("{} \\".format(sys.argv[0]))
        for option, value in vars(args).items():
            option = option.replace("_", "-")
            print("    --{} {} \\".format(option, value))
        print()

        # Checking the input files
        input_files = check_input_files(args.input)

        # Getting the tool's configuration and setting it
        tool_config = read_config_file(args.tool_config)
        Tool.set_tool_configuration(tool_config)

        # Do we run using DRMAA?
        if args.use_drmaa:
            Tool.do_not_run_locally()

        # Getting the pipeline steps
        what_to_run = get_pipeline_steps(args.pipeline_config)

        # The first step of the pipeline
        @originate(input_files)
        def start(o_files):
            pass

        # Dynamically creating the pipeline
        job_order = []
        in_job = start
        last_suffix = ""
        curr_formatter = None
        curr_output = None
        for job_index, (job, job_options) in enumerate(what_to_run):
            # Getting the input and output file type
            input_type = job.get_input_type()
            output_type = job.get_output_type()

            # The output directory
            output_dir = os.path.join("output",
                                    "{:02d}_{}".format(job_index + 1,
                                                        job.get_tool_name()))
            # Creating the output directory
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)

            # What will be in the formatter
            curr_formatter = [r".+/(?P<SAMPLE>[a-zA-Z0-9_\-]+){}".format(i)
                                                            for i in input_type]
            curr_output = [os.path.join(output_dir, ("{SAMPLE[" + str(i) +"]}" +
                                                    last_suffix + suffix))
                                        for i, suffix in enumerate(output_type)]
            formatter_func = formatter

            # What if we need to merge all inputs?
            if job.need_to_merge_all_inputs():
                curr_formatter = [r".+/[a-zA-Z0-9_\-]+{}".format(i)
                                                            for i in input_type]
                curr_output = [os.path.join(output_dir, ("all_samples" +
                                                        last_suffix + suffix))
                                                      for suffix in output_type]
                formatter_func = regex

            # Checking if there is only one output
            if len(curr_output) == 1:
                curr_output = curr_output[0]

            # Getting the current Ruffus' decorator
            curr_decorator = None
            if (len(input_type) > len(output_type)) or job.need_to_merge_all_inputs():
                curr_decorator = collate
            elif len(input_type) == len(output_type):
                curr_decorator = transform
            else:
                m = "cannot choose a good Ruffus decorator"
                raise ProgramError(m)

            # The name of the function
            func_name = "step{:02d}_{}".format(job_index + 1,
                                               job.get_tool_name())

            # Dynamically creating the pipeline
            @curr_decorator(in_job, formatter_func(*curr_formatter),
                            curr_output, "{SAMPLE[0]}", job, len(input_type),
                            len(output_type), output_dir, job_options)
            @rename_func(func_name)
            def curr_step(i_files, o_files, sample_id, job, nb_in, nb_out,
                          out_dir, options):
                print("\n###########################")
                print(job.get_tool_name())
                # If we need to merge all inputs
                if job.need_to_merge_all_inputs():
                    i_files = list(i_files)
                    sample_id = "all_samples"

                # The i_files variable is usually a tuple of lists
                if isinstance(i_files, tuple):
                    i_files = i_files[0]

                # We want to work on a copy of the options
                curr_options = copy(options)

                # Adding the input to the tool option
                if job.need_to_merge_all_inputs():
                    curr_options["inputs"] = i_files
                elif nb_in == 1:
                    curr_options["input"] = i_files
                else:
                    for i in range(nb_in):
                        curr_options["input{}".format(i + 1)] = i_files[i]

                # Adding the output files
                if nb_out == 1:
                    curr_options["output"] = o_files
                else:
                    for i in range(nb_out):
                        curr_options["output{}".format(i + 1)] = o_files[i]

                # Adding the prefix and sample id
                if "prefix" not in curr_options:
                    curr_options["prefix"] = os.path.join(out_dir, sample_id)
                if "sample_id" not in curr_options:
                    curr_options["sample_id"] = sample_id
                print(curr_options["sample_id"])

                # Running the task
                job.execute(curr_options, out_dir=out_dir)

            # Setting the attribute for the new function so that it can be
            # pickled
            setattr(__main__, func_name, curr_step)

            # Updating the in_job and the last suffix only if the tool produces
            # usable data
            if job.produce_usable_data():
                in_job = curr_step
                last_suffix += ".{}".format(job.get_suffix())

            # Adding the current job to the pipeline
            job_order.append(curr_step)

        # Printing the pipeline
        print("Running the pipeline...")
        pipeline_printout_graph("flowchart.svg", "svg", job_order)
        pipeline_run(job_order, verbose=0, multiprocess=args.nb_process,
                     checksum_level=1)

    except KeyboardInterrupt:
        print("Cancelled by user", sys.stderr)
        sys.exit(0)
    except ProgramError as e:
        parser.error(e.message)
