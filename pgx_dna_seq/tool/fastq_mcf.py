__all__ = ["ClipTrim"]

from pgx_dna_seq.tool import GenericTool


class FastQMCF(GenericTool):

    # The version of the tool
    _version = "1.1.2-537"

    # The executable
    _exec = "fastq-mcf"

    def __init__(self):
        """Initialize a FastQMCF instance."""
        pass


class ClipTrim(FastQMCF):

    # The name of the tool
    _tool_name = "ClipTrim"

    # The options
    _command = ("{adapters} {input1} {input2} -o {output1} -o {output2} "
                "{other_opt}")

    # The STDOUT and STDERR
    _stdout = "{prefix}.out"
    _stderr = "{prefix}.err"

    # The description of the required options
    _required_options = {"input1":    GenericTool.INPUT,
                         "input2":    GenericTool.INPUT,
                         "output1":   GenericTool.OUTPUT,
                         "output2":   GenericTool.OUTPUT,
                         "adapters":  GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "prefix":    GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "ct"

    # The input and output type
    _input_type = (r"_R1\.(\S+\.)?fastq(\.gz)?$",
                   r"_R2\.(\S+\.)?fastq(\.gz)?$")
    _output_type = ("_R1.{}.fastq.gz".format(_suffix),
                    "_R2.{}.fastq.gz".format(_suffix))

    def __init__(self):
        """Initialize a ClipTrim instance."""
        pass
