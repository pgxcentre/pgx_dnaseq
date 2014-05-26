__all__ = ["CoverageGraph", "CoverageGraph_Multi"]

import re

from pgx_dna_seq.tool import GenericTool
from pgx_dna_seq.tool.samtools import IndexBam


class PGx_CoverageGraph(GenericTool):

    # The version of the tool
    _version = "0.1"

    # The executable
    _exec = "coverage_graph.py"

    def __init__(self):
        """Initialize a PGx_CoverageGraph instance."""
        pass


class CoverageGraph(PGx_CoverageGraph):

    # The name of the tool
    _tool_name = "CoverageGraph"

    # The options
    _command = "{other_opt} --out {out_prefix} --bam {input}"

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":      GenericTool.INPUT,
                         "output":     GenericTool.OUTPUT,
                         "out_prefix": GenericTool.OUTPUT,
                         "other_opt":  GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
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
        super(CoverageGraph, self).execute(options, out_dir)


class CoverageGraph_Multi(PGx_CoverageGraph):

    # The name of the tool
    _tool_name = "CoverageGraph_Multi"

    # The options
    _command = "{other_opt} --out {out_prefix} --bam {inputs}"

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"inputs":     GenericTool.INPUTS,
                         "output":     GenericTool.OUTPUT,
                         "out_prefix": GenericTool.OUTPUT,
                         "other_opt":  GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
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
        super(CoverageGraph_Multi, self).execute(options, out_dir)
