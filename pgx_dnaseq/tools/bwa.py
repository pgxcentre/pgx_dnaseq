
# This file is part of pgx_dnaseq
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os

from . import GenericTool
from .. import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["ALN", "SAMPE", "MEM"]


class BWA(GenericTool):

    # The version of the TOOL
    _version = "0.7.12"

    # The executable
    _exec = "bwa"

    def __init__(self):
        """Initialize a BWA instance."""
        pass


class ALN(BWA):

    # The name of the tool
    _tool_name = "ALN"

    # The options
    _command = "aln {reference} {other_aln_opt} {input}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"reference":     GenericTool.INPUT,
                         "other_aln_opt": GenericTool.OPTIONAL,
                         "input":         GenericTool.INPUT,
                         "output":        GenericTool.OUTPUT}

    def __init__(self):
        """Initialize a ALN instance."""
        pass


class MEM(BWA):

    # The name of the tool
    _tool_name = "MEM"

    # The options
    _command = "mem {reference} {input1} {input2} {other_opt} "

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "input1":    GenericTool.INPUT,
                         "input2":    GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "mem"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = (".{}.sam".format(_suffix),)

    def __init__(self):
        """Initialize a MEM instance."""
        pass


class SAMPE(BWA):

    # The name of the tool
    _tool_name = "SAMPE"

    # The options
    _command = "sampe {reference} {sai1} {sai2} {input1} {input2}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"reference": GenericTool.INPUT,
                         "sai1":      GenericTool.INPUT,
                         "sai2":      GenericTool.INPUT,
                         "input1":    GenericTool.INPUT,
                         "input2":    GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "sampe"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = (".{}.sam".format(_suffix),)

    def __init__(self):
        """Initialize a SAMPE instance."""
        pass

    def execute(self, options, out_dir=None):
        """Execute ALN and SAMPE."""
        # The ALN options for the first file
        aln_options = {}
        try:
            # The reference
            aln_options["reference"] = options["reference"]

            # The input file
            aln_options["input"] = options["input1"]

            # the output SAI
            input_basename = os.path.basename(aln_options["input"])
            output_sai = "{}.sai".format(input_basename)
            output_sai = os.path.join(os.path.dirname(options["output"]),
                                      output_sai)
            aln_options["output"] = output_sai
            options["sai1"] = aln_options["output"]
            if "other_aln_opt" in options:
                aln_options["other_aln_opt"] = options["other_aln_opt"]
        except KeyError as e:
            m = "{}: missing option {}".format(self.__class__.__name__, e)
            raise ProgramError(m)

        # Execution for the first file
        ALN().execute(aln_options, out_dir)

        # The ALN options for the second file
        try:
            # The input file
            aln_options["input"] = options["input2"]

            # The output SAI
            input_basename = os.path.basename(aln_options["input"])
            output_sai = "{}.sai".format(input_basename)
            output_sai = os.path.join(os.path.dirname(options["output"]),
                                      output_sai)
            aln_options["output"] = output_sai
            options["sai2"] = aln_options["output"]
        except KeyError as e:
            m = "{}: missing option {}".format(self.__class__.__name__, e)
            raise ProgramError(m)

        # Executing for the second file
        ALN().execute(aln_options, out_dir)

        # Executing SAMPE
        super().execute(options, out_dir)
