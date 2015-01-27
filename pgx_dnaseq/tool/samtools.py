
# This file is part of pgx_dnaseq
#
# This work is licensed under The MIT License (MIT). To view a copy of this
# license, visit http://opensource.org/licenses/MIT


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "MIT"


import re

from . import GenericTool
from .. import ProgramError


__all__ = ["Sam2Bam", "IndexBam", "KeepMapped", "FlagStat", "MPILEUP",
           "MPILEUP_Multi"]


class Samtools(GenericTool):

    # The version of the tool
    _version = "1.1"

    # The executable
    _exec = "samtools"

    def __init__(self):
        """Initialize a Samtools instance."""
        pass


class Sam2Bam(Samtools):

    # The name of the tool
    _tool_name = "Sam2Bam"

    # The options
    _command = "view -h -b -S {input}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":  GenericTool.INPUT,
                         "output": GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "sam2bam"

    # The input and output type
    _input_type = (r"\.(\S+\.)?sam$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
        """Initialize a Sam2Bam instance."""
        pass


class IndexBam(Samtools):

    # The name of the tool
    _tool_name = "IndexBam"

    # The options
    _command = "index {input}"

    # The STDOUT and STDERR
    _stdout = "{input}.bai.out"
    _stderr = "{input}.bai.err"

    # The description of the required options
    _required_options = {"input": GenericTool.INPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = None

    # The input and output type
    _input_type = (r"\.(\S+\.)?sam$", )
    _output_type = (".bai", )

    def __init__(self):
        """Initialize a IndexBam instance."""
        pass


class KeepMapped(Samtools):

    # The name of the tool
    _tool_name = "KeepMapped"

    # The options
    _command = "view -u -h {mapped_opt} {input}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":      GenericTool.INPUT,
                         "output":     GenericTool.OUTPUT,
                         "mapped_opt": GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "mapped"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
        """Initialize a KeepMapped instance."""
        pass

    def execute(self, options, out_dir=None):
        """Extract mapped reads (keeping the unmapped ones)."""
        # First, we extract the mapped reads
        options["mapped_opt"] = "-F 4"
        super(KeepMapped, self).execute(options, out_dir)

        # Then, we extract the unmapped reads
        options["mapped_opt"] = "-f 4"
        options["output"] = re.sub("{}$".format(self._output_type[0]),
                                   ".unmapped.bam", options["output"])
        super(KeepMapped, self).execute(options, out_dir)


class FlagStat(Samtools):

    # The name of the tool
    _tool_name = "FlagStat"

    # The options
    _command = "flagstat {input}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":  GenericTool.INPUT,
                         "output": GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "flagstat"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}".format(_suffix), )

    # This tool does not produce usable data...
    _produce_data = False

    def __init__(self):
        """Initialize a FlagStat instance."""
        pass


class MPILEUP(Samtools):

    # The name of the tool
    _tool_name = "MPILEUP"

    # The options
    _command = "mpileup {other_opt} -f {reference} {input}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input": GenericTool.INPUT,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "output": GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "mpileup"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}".format(_suffix), )

    def __init__(self):
        """Initialize a MPILEUP instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes a BAM and create the MPILEUP file."""
        # First we index the input file
        if "input" not in options:
            m = "{}: no input file".format(self.__class__.__name__)
            raise ProgramError(m)
        IndexBam().execute({"input": options["input"]}, out_dir)

        # Then we create the MPILEUP file
        super(MPILEUP, self).execute(options, out_dir)


class MPILEUP_Multi(Samtools):

    # The name of the tool
    _tool_name = "MPILEUP_Multi"

    # The options
    _command = "mpileup {other_opt} -f {reference} {inputs}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"inputs": GenericTool.INPUTS,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "output": GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "mpileup"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}".format(_suffix), )

    # This tool needs multiple input
    _merge_all_inputs = True

    def __init__(self):
        """Initialize a MPILEUP_Multi instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes a BAM and create the MPILEUP file."""
        # First we index all the input files
        if "inputs" not in options:
            m = "{}: no input files".format(self.__class__.__name__)
            raise ProgramError(m)
        for filename in options["inputs"]:
            IndexBam().execute({"input": filename}, out_dir)

        # Then we create the MPILEUP file from multiple inputs
        super(MPILEUP_Multi, self).execute(options, out_dir)
