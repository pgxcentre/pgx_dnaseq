
# This file is part of pgx_dnaseq
#
# This work is licensed under The MIT License (MIT). To view a copy of this
# license, visit http://opensource.org/licenses/MIT


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "MIT"


import re
from glob import glob

from . import GenericTool


__all__ = ["ClipTrim"]


class FastQMCF(GenericTool):

    # The version of the tool
    _version = "1.1.2-806"

    # The executable
    _exec = "fastq-mcf"

    def __init__(self):
        """Initialize a FastQMCF instance."""
        pass


class ClipTrim(FastQMCF):

    # The name of the tool
    _tool_name = "ClipTrim"

    # The options
    _command = ("{adapters} {input1} {input2} -o {output1} -o {output2} "
                "{other_opt}")

    # The STDOUT and STDERR
    _stdout = "{prefix}.out"
    _stderr = "{prefix}.err"

    # The description of the required options
    _required_options = {"input1":    GenericTool.INPUT,
                         "input2":    GenericTool.INPUT,
                         "output1":   GenericTool.OUTPUT,
                         "output2":   GenericTool.OUTPUT,
                         "adapters":  GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "prefix":    GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "ct"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = ("_R1.{}.fastq.gz".format(_suffix),
                    "_R2.{}.fastq.gz".format(_suffix))

    def __init__(self):
        """Initialize a ClipTrim instance."""
        pass

    def read_report(self, prefix):
        """Reads a ClipTrim report file."""
        # Getting the report file name
        filename = glob("{}.out".format(prefix))
        assert len(filename) == 1
        filename = filename[0]

        # Reading the content
        content = None
        with open(filename, "r") as i_file:
            content = i_file.read()

        # Getting the total number of reads
        nb_reads = int(re.search(r"Total reads: (\d+)", content).group(1)) * 2

        # Getting the number of reads that were too short
        nb_shorts = re.search(r"Too short after clip: (\d+)", content)
        nb_shorts = int(nb_shorts.group(1)) * 2

        # The number of trimmed reads
        trim_r1 = re.search(r"Trimmed (\d+) reads .*_R1\.fastq\.gz", content)
        trim_r1 = int(trim_r1.group(1))
        trim_r2 = re.search(r"Trimmed (\d+) reads .*_R2\.fastq\.gz", content)
        trim_r2 = int(trim_r2.group(1))

        # Saving the results
        result = {
            "total_reads_before_trim":   nb_reads,
            "nb_short_reads_after_trim": nb_shorts,
            "nb_trimmed_r1":             trim_r1,
            "nb_trimmed_r2":             trim_r2,
        }

        return result

