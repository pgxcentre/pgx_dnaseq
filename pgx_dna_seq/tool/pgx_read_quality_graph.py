__all__ = ["ReadQualityGraph", ]

from glob import glob
from shutil import copyfile

from pgx_dna_seq.tool import GenericTool


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
                "--log {output}.log")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input1":     GenericTool.INPUT,
                         "input2":     GenericTool.INPUT,
                         "output":     GenericTool.OUTPUT,
                         "other_opt":  GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
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

