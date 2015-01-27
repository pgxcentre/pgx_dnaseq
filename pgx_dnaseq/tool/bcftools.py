
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


__all__ = ["BcftoolsVariantCaller"]


class Bcftools(GenericTool):

    # The version of the tool
    _version = "1.1"

    # The executable
    _exec = "bcftools"

    def __init__(self):
        """Initialize a Bcftools instance."""
        pass


class BcftoolsVariantCaller(Bcftools):

    # The name of the tool
    _tool_name = "BcftoolsVariantCaller"

    # The options
    _command = "view -vcg {input}"

    # The STDOUT and STDERR
    _stdout = "{output}"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":  GenericTool.INPUT,
                         "output": GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "bcftools"

    # The input and output type
    _input_type = (r"\.(\S+\.)?mpileup$", )
    _output_type = (".{}.vcf".format(_suffix), )

    def __init__(self):
        """Initialize a BcftoolsVariantCaller instance."""
        pass


