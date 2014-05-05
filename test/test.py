#!/usr/bin/env python

import os
import sys
import __main__

from ruffus import *
from pgx_dna_seq.tool.bwa import SAMPE
from pgx_dna_seq.tool.fastq_mcf import ClipTrim
from pgx_dna_seq.tool.fastqc import FastQC_FastQ
from pgx_dna_seq.tool.gatk import BaseRecalibrator, IndelRealigner
from pgx_dna_seq.tool.samtools import Sam2Bam, FlagStat, KeepMapped
from pgx_dna_seq.tool.picard_tools import AddRG, SortSam, MarkDuplicates


# Creating a decorator that changes a function name
def rename_func(new_name):
    """Renames a function."""
    def decorator(func):
        func.__name__ = new_name
        return func
    return decorator

# The input files
starting_files = [["data/NA12877_R1.fastq.gz", "data/NA12877_R2.fastq.gz"],
                  ["data/NA12878_R1.fastq.gz", "data/NA12878_R2.fastq.gz"]]

# The first step
@originate(starting_files)
def start(o_files):
    print(o_files)
    pass

# What to run
what_to_run = [FastQC_FastQ(), ClipTrim(), FastQC_FastQ(), SAMPE(), Sam2Bam(),
               FlagStat(), KeepMapped(), SortSam(), AddRG(), IndelRealigner(),
               BaseRecalibrator(), MarkDuplicates()]

# The options
options = [{},                                      # FastQC

           {"adapters":  "data/adapters.fa",        # Clip/Trim
            "other_opt": "-S -q 20"},

           {},                                      # FastQC

           {"reference":     "reference/hg19.fasta",# SAMPE
            "other_aln_opt": "-t 4"},

           {},                                      # Sam 2 Bam

           {},                                      # Flagstat

           {},                                      # Keep mapped

           {"java_memory":    "6g",                 # Sort
            "java_other_opt": "-Djava.io.tmpdir=TMP",
            "sort_order":     "coordinate",
            "other_options":  "TMP_DIR=TMP"},

           {"rgpl":      "ILLUMINA",                # Add RG
            "rgds":      "SureSelectExomeV4_HiSeq2000",
            "rgpu":      "Candid",
            "rgcn":      "Universite_de_Montreal_Centre_de_Pharmacogenomique_Beaulieu-Saucier",
            "rglb":      "run_library"},

           {"reference":   "reference/hg19.fasta",            # Realign
            "dt":          "NONE",
            "java_memory": "6g"},

           {"java_memory":       "14g",             # Recal
            "reference":         "reference/hg19.fasta",
            "dt":                "NONE",
            "dbSNP_known_sites": "reference/dbSNP_138.GRCh37_p10.vcf.gz"},

           {"java_memory": "10g"}]                  # Mark duplicate

# The options for DRMAA
drmaa_options = {"ALN":                    {"walltime": "00:15:00",
                                            "nb_node":  1,
                                            "nb_proc":  4},

                 "SAMPE":                  {"walltime": "00:15:00",
                                            "nb_node":  1,
                                            "nb_proc":  4},

                 "ClipTrim":               {"walltime": "00:06:00",
                                            "nb_node":  1,
                                            "nb_proc":  1},

                 "FastQC_FastQ":           {"walltime": "00:06:00",
                                            "nb_node":  1,
                                            "nb_proc":  1},

                 "RealignerTargetCreator": {"walltime": "00:60:00",
                                            "nb_node":  1,
                                            "nb_proc":  8},

                 "IndelRealigner":         {"walltime": "00:60:00",
                                            "nb_node":  1,
                                            "nb_proc":  8},

                 "PrintReads":             {"walltime": "00:60:00",
                                            "nb_node":  1,
                                            "nb_proc":  4},

                 "BaseRecalibrator":       {"walltime": "00:60:00",
                                            "nb_node":  1,
                                            "nb_proc":  4},

                 "SortSam":                {"walltime": "00:15:00",
                                            "nb_node":  1,
                                            "nb_proc":  4},

                 "AddRG":                  {"walltime": "00:06:00",
                                            "nb_node":  1,
                                            "nb_proc":  1},

                 "MarkDuplicates":         {"walltime": "00:15:00",
                                            "nb_node":  1,
                                            "nb_proc":  8},

                 "Sam2Bam":                {"walltime": "00:06:00",
                                            "nb_node":  1,
                                            "nb_proc":  1},

                 "IndexBam":               {"walltime": "00:06:00",
                                            "nb_node":  1,
                                            "nb_proc":  1},

                 "KeepMapped":             {"walltime": "00:06:00",
                                            "nb_node":  1,
                                            "nb_proc":  1},

                 "FlagStat":               {"walltime": "00:06:00",
                                            "nb_node":  1,
                                            "nb_proc":  1},}

# The job order
job_order = []

# Dynamically creating the pipeline
in_job = start
last_suffix = ""
curr_formatter = None
curr_output = None
for curr_job, job in enumerate(what_to_run):
    # Getting the input and output file type
    input_type = job.get_input_type()
    output_type = job.get_output_type()

    # The output directory
    output_dir = os.path.join("output",
                              "{:02d}_{}".format(curr_job + 1,
                                                 job.get_tool_name()))
    # Creating the output directory
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # What will be in the formatter
    curr_formatter = [r".+/(?P<SAMPLE>\w+){}".format(i) for i in input_type]
    curr_output = [os.path.join(output_dir, ("{SAMPLE[" + str(i) +"]}" +
                                             last_suffix + suffix))
                                for i, suffix in enumerate(output_type)]

    print("\n".join(curr_formatter))
    print("\n".join(curr_output))
    print()

    # Checking if there is only one output
    if len(curr_output) == 1:
        curr_output = curr_output[0]

    # Getting the current Ruffus' decorator
    curr_decorator = None
    if len(input_type) == len(output_type):
        curr_decorator = transform
    elif len(input_type) > len(output_type):
        curr_decorator = collate

    # The name of the function
    func_name = "step{:02d}_{}".format(curr_job + 1, job.get_tool_name())

    # Dynamically creating the pipeline
    @curr_decorator(in_job, formatter(*curr_formatter), curr_output,
                    "{SAMPLE[0]}", job, len(input_type), len(output_type),
                    output_dir, options[curr_job], drmaa_options)
    @rename_func(func_name)
    def curr_step(i_files, o_files, sample_id, job, nb_in, nb_out, out_dir,
                  options, drmaa_options):
        print("\n###########################")
        print(job.get_tool_name())
        # The i_files variable is a tuple of lists
        if isinstance(i_files, tuple):
            i_files = i_files[0]

        print(i_files)
        print(o_files)
        print()

        # Adding the input to the tool option
        if nb_in == 1:
            options["input"] = i_files
        else:
            for i in range(nb_in):
                options["input{}".format(i + 1)] = i_files[i]

        # Adding the output files
        if nb_out == 1:
            options["output"] = o_files
        else:
            for i in range(nb_out):
                options["output{}".format(i + 1)] = o_files[i]

        # Adding the prefix and sample id
        if "prefix" not in options:
            options["prefix"] = os.path.join(out_dir, sample_id)
        if "sample_id" not in options:
            options["sample_id"] = sample_id

        # Running the task
        job.execute(options, drmaa_options=drmaa_options, out_dir=out_dir,
                    locally=True)

    # Setting the attribute for the new function so that it can be pickled
    setattr(__main__, func_name, curr_step)

    # Updating the in_job and the last suffix only if the tool produces usable
    # data
    if job.produce_usable_data():
        in_job = curr_step
        last_suffix += ".{}".format(job.get_suffix())
    else:
        in_job = output_from([in_job, curr_step])


print("Runing the pipeline...")
pipeline_printout_graph("flowchart.svg", "svg", job_order)
pipeline_run(verbose=0, multiprocess=1, checksum_level=1)
pipeline_printout_graph("flowchart_after.svg", "svg", job_order)
