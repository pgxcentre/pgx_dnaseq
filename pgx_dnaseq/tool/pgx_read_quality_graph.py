
# This file is part of pgx_dnaseq
#
# This work is licensed under The MIT License (MIT). To view a copy of this
# license, visit http://opensource.org/licenses/MIT


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "MIT"


from glob import glob
from shutil import copyfile

from . import GenericTool


__all__ = ["ReadQualityGraph", ]


class PGx_ReadQualityGraph(GenericTool):

    # The version of the tool
    _version = "0.1"

    # The executable
    _exec = "read_quality_graph.py"

    def __init__(self):
        """Initialize a PGx_ReadQualityGraph instance."""
        pass


class ReadQualityGraph(PGx_ReadQualityGraph):

    # The name of the tool
    _tool_name = "ReadQualityGraph"

    # The options
    _command = ("{other_opt} --input {input1} {input2} --output {output} "
                "--title-prefix {sample_id} --log {output}.log")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input1":    GenericTool.INPUT,
                         "input2":    GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "sample_id": GenericTool.REQUIREMENT,
                         "other_opt": GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output
    # file
    _suffix = "read_quality_graph"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = (".{}.pdf".format(_suffix),)

    # This tool does not produce usable data...
    _produce_data = False

    def __init__(self):
        """Initialize a ReadQualityGraph instance."""
        pass

    def read_report(self, prefix):
        """Reads a ReadQualityGraph report file."""
        # Getting the report file name
        filename = glob("{}*.{}.pdf".format(prefix, self._suffix))
        assert len(filename) == 1
        filename = filename[0]

        # Copying the file
        new_filename = "{}_{}.pdf".format(prefix, self.get_tool_name())
        copyfile(filename, new_filename)

        # Saving the results
        result = {
            "read_qual_figname": new_filename,
        }

        return result
