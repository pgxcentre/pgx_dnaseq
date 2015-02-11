
# This file is part of pgx_dnaseq
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import re
from glob import glob
from shutil import copyfile

from . import GenericTool
from .samtools import IndexBam


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["CoverageGraph", "CoverageGraph_Multi"]


class PGx_CoverageGraph(GenericTool):

    # The version of the tool
    _version = "0.3"

    # The executable
    _exec = "coverage_graph.py"

    def __init__(self):
        """Initialize a PGx_CoverageGraph instance."""
        pass


class CoverageGraph(PGx_CoverageGraph):

    # The name of the tool
    _tool_name = "CoverageGraph"

    # The options
    _command = "{other_opt} --out {out_prefix} --bed {targets} --bam {input}"

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":      GenericTool.INPUT,
                         "output":     GenericTool.OUTPUT,
                         "out_prefix": GenericTool.OUTPUT,
                         "targets":    GenericTool.INPUT,
                         "other_opt":  GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "coverage_graph"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.png".format(_suffix), )

    # This tool does not produce usable data...
    _produce_data = False

    def __init__(self):
        """Initialize a CoverageGraph instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes all the BAM files and plots the coverage."""
        # First we index all the input files
        if "input" not in options:
            m = "{}: no input file".format(self.__class__.__name__)
            raise ProgramError(m)
        IndexBam().execute({"input": options["input"]}, out_dir)

        # Searching for the out prefix
        out_prefix = re.sub("\.png$", "", options["output"])
        options["out_prefix"] = out_prefix

        # Executing the software
        super().execute(options, out_dir)

    def read_report(self, prefix):
        """Reads a CoverageGraph report file."""
        # Getting the report file name
        filename = glob("{}*.{}.png".format(prefix, self._suffix))
        assert len(filename) == 1
        filename = filename[0]

        # Copying the file
        new_filename = "{}_{}.png".format(prefix, self.get_tool_name())
        copyfile(filename, new_filename)

        # Saving the results
        result = {
            "coverage_figname": new_filename,
        }

        return result


class CoverageGraph_Multi(PGx_CoverageGraph):

    # The name of the tool
    _tool_name = "CoverageGraph_Multi"

    # The options
    _command = "{other_opt} --out {out_prefix} --bed {targets} --bam {inputs}"

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"inputs":     GenericTool.INPUTS,
                         "output":     GenericTool.OUTPUT,
                         "out_prefix": GenericTool.OUTPUT,
                         "targets":    GenericTool.INPUT,
                         "other_opt":  GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "coverage_graph"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.png".format(_suffix), )

    # This tool needs multiple input
    _merge_all_inputs = True

    # This tool does not produce usable data...
    _produce_data = False

    def __init__(self):
        """Initialize a CoverageGraph_Multi instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes all the BAM files and plots the coverage."""
        # First we index all the input files
        if "inputs" not in options:
            m = "{}: no input files".format(self.__class__.__name__)
            raise ProgramError(m)
        for filename in options["inputs"]:
            IndexBam().execute({"input": filename}, out_dir)

        # Searching for the out prefix
        out_prefix = re.sub("\.png$", "", options["output"])
        options["out_prefix"] = out_prefix

        # Executing the software
        super().execute(options, out_dir)

    def read_report(self, prefix):
        """Reads a CoverageGraph report file."""
        # Getting the report file name
        filename = glob("{}*.{}.png".format(prefix, self._suffix))
        assert len(filename) == 1
        filename = filename[0]

        # Copying the file
        new_filename = "{}_{}.png".format(prefix, self.get_tool_name())
        copyfile(filename, new_filename)

        # Saving the results
        result = {
            "coverage_multi_figname": new_filename,
        }

        return result
