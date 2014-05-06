__all__ = ["Sam2Bam", "IndexBam", "KeepMapped", "FlagStat"]

import re

from pgx_dna_seq.tool import GenericTool


class Samtools(GenericTool):

    # The version of the tool
    _version = "0.1.19"

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

    # The suffix that will be added just before the extension of the output file
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

    # The suffix that will be added just before the extension of the output file
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

    # The suffix that will be added just before the extension of the output file
    _suffix = "mapped"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
        """Initialize a KeepMapped instance."""
        pass

    def execute(self, options, out_dir=None, locally=True):
        """Extract mapped reads (keeping the unmapped ones)."""
        # First, we extract the mapped reads
        options["mapped_opt"] = "-F 4"
        super(KeepMapped, self).execute(options, out_dir, locally)

        # Then, we extract the unmapped reads
        options["mapped_opt"] = "-f 4"
        options["output"] = re.sub("{}$".format(self._output_type[0]),
                                   ".unmapped.bam", options["output"])
        super(KeepMapped, self).execute(options, out_dir, locally)


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

    # The suffix that will be added just before the extension of the output file
    _suffix = "flagstat"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}".format(_suffix), )

    # This tool does not produce usable data...
    _produce_data = False

    def __init__(self):
        """Initialize a FlagStat instance."""
        pass
