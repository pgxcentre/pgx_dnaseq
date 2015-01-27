
# This file is part of pgx_dnaseq
#
# This work is licensed under The MIT License (MIT). To view a copy of this
# license, visit http://opensource.org/licenses/MIT


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "MIT"


import os
import re
import shutil

from . import GenericTool
from .. import ProgramError


__all__ = ["FastQC_FastQ"]


class FastQC(GenericTool):

    # The version of the tool
    _version = "0.11.2"

    # The executable
    _exec = "fastqc"

    def __init__(self):
        """Initialize a FastQC instance."""
        pass


class FastQC_FastQ(FastQC):

    # The name of the tool
    _tool_name = "FastQC_FastQ"

    # The options
    _command = "{input}"

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":  GenericTool.INPUT,
                         "output": GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "fastqc"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = ("_R1_{}.zip".format(_suffix),
                    "_R2_{}.zip".format(_suffix))

    # This tool does not produce usable data...
    _produce_data = False

    def __init__(self):
        """Initialize a FastQC_FastQ instance."""
        pass

    def execute(self, options, out_dir=None):
        """Execute FastQC on both input files."""
        # First, we run on the first file
        try:
            options["input"] = options["input1"]
            options["output"] = options["output1"]
        except KeyError as e:
            m = "{}: missing option {}".format(self.__class__.__name__, e)
            raise ProgramError(m)
        super(FastQC_FastQ, self).execute(options, out_dir)

        # We need to move the files... FastQC renames the original file so that
        # the directory is the same, but ".fastq(.gz)?" is replaced to
        # "_fastqc" (for the directory) and "_fastqc.zip" for the file.
        dir_to_move = re.sub(r"\.fastq(\.gz)$", "_fastqc", options["input"])
        destination = re.sub(r"\.zip$", "", options["output1"])
        if os.path.isdir(dir_to_move):
            shutil.move(dir_to_move, destination)
        file_to_move = "{}.zip".format(dir_to_move)
        if os.path.isfile(file_to_move):
            shutil.move(file_to_move, options["output1"])

        # Then, we run on the second file
        try:
            options["input"] = options["input2"]
            options["output"] = options["output2"]
        except KeyError as e:
            m = "{}: missing option {}".format(self.__class__.__name__, e)
            raise ProgramError(m)
        super(FastQC_FastQ, self).execute(options, out_dir)

        # We need to move the files... FastQC renames the original file so that
        # the directory is the same, but ".fastq(.gz)?" is replaced to
        # "_fastqc" (for the directory) and "_fastqc.zip" for the file.
        dir_to_move = re.sub(r"\.fastq(\.gz)$", "_fastqc", options["input"])
        destination = re.sub(r"\.zip$", "", options["output2"])
        if os.path.isdir(dir_to_move):
            shutil.move(dir_to_move, destination)
        file_to_move = "{}.zip".format(dir_to_move)
        if os.path.isfile(file_to_move):
            shutil.move(file_to_move, options["output2"])
