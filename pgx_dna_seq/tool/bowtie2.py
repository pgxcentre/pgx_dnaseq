__all__ = ["Bowtie2_align"]

from pgx_dna_seq.tool import GenericTool


class Bowtie2(GenericTool):

    # The version of the tool
    _version = "2.2.2"

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

    # The suffix that will be added just before the extension of the output file
    _suffix = "bowtie2"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = (".{}.sam".format(_suffix),)

    def __init__(self):
        """Initialize a Bowtie2_align instance."""
        pass
