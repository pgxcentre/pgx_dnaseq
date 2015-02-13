
# This file is part of pgx_dnaseq
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import re

from . import GenericTool
from .. import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["VcfConcat", "VcfSort"]


class Vcftools(GenericTool):

    # The version of the tool
    _version = "0.1.12b"

    # The executable
    _exec = "vcftools"

    def __init__(self):
        """Initialize a Vcftools instance."""
        pass


class VcfConcat(Vcftools):

    # The name of the tool
    _tool_name = "VcfConcat"

    # The executable is different
    _exec = "vcf-concat"

    # The options
    _command = "{other_opt} {inputs}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"inputs":    GenericTool.INPUTS,
                         "other_opt": GenericTool.OPTIONAL,
                         "output":    GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "vcf_concat"

    # The input and output type
    _input_type = (r"\.(\S+\.)?vcf$", )
    _output_type = (".{}.vcf".format(_suffix), )

    def __init__(self):
        """Initialize a VcfConcat instance."""
        pass


class VcfSort(Vcftools):

    # The name of th etool
    _tool_name = "VcfSort"

    # The executable is different
    _exec = "vcf-sort"

    # The options
    _command = "{other_opt} {input}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "output":    GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "vcf_sort"

    # The input and output type
    _input_type = (r"\.(\S+\.)?vcf$", )
    _output_type = (".{}.vcf".format(_suffix), )

    def __init__(self):
        """Initialize a VcfSort instance."""
        pass
