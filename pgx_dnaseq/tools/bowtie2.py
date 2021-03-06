
# This file is part of pgx_dnaseq
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


from . import GenericTool


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["Bowtie2_align"]


class Bowtie2(GenericTool):

    # The version of the tool
    _version = "2.2.4"

    # The executable
    _exec = "bowtie2"

    def __init__(self):
        """Initialize a Bowtie2  instance."""
        pass


class Bowtie2_align(Bowtie2):

    # The name of the tool
    _tool_name = "Bowtie2_align"

    # The options
    _command = "{reference} -1 {input1} -2 {input2} -S {output} {other_opt}"

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input1":    GenericTool.INPUT,
                         "input2":    GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.REQUIREMENT,
                         "other_opt": GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "bowtie2"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = (".{}.sam".format(_suffix),)

    def __init__(self):
        """Initialize a Bowtie2_align instance."""
        pass
