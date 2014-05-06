__all__ = ["SortSam", "AddRG", "MarkDuplicates"]

import re

from pgx_dna_seq.tool.java import JAR
from pgx_dna_seq.tool import GenericTool


class PicardTools(JAR):

    # The version of the tool
    _version = "1.113"

    def __init__(self):
        """Initialize a PicardTools instance."""
        pass


class SortSam(PicardTools):

    # The name of the tool
    _tool_name = "SortSam"

    # The jar location
    _jar = "/opt/picard-tools-1.113/SortSam.jar"

    # The options
    _command = ("INPUT={input} OUTPUT={output} SORT_ORDER={sort_order} "
                "{other_options}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":         GenericTool.INPUT,
                         "output":        GenericTool.OUTPUT,
                         "sort_order":    GenericTool.REQUIREMENT,
                         "other_options": GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
    _suffix = "sorted"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
       """Initialize a SortSam instance."""
       pass


class AddRG(PicardTools):

    # The name of the tool
    _tool_name = "AddRG"

    # The jar location
    _jar = "/opt/picard-tools-1.113/AddOrReplaceReadGroups.jar"

    # The options
    _command = ("I={input} O={output} RGDS={rgds} RGPL={rgpl} RGPU={rgpu} "
                "RGSM={sample_id} RGCN={rgcn} RGLB={rglb}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":       GenericTool.INPUT,
                         "output":      GenericTool.OUTPUT,
                         "rgds":        GenericTool.REQUIREMENT,
                         "rgpl":        GenericTool.REQUIREMENT,
                         "rgpu":        GenericTool.REQUIREMENT,
                         "sample_id":   GenericTool.REQUIREMENT,
                         "rgcn":        GenericTool.REQUIREMENT,
                         "rglb":        GenericTool.REQUIREMENT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "rg"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
        """Initialize a AddRG instance."""
        pass


class MarkDuplicates(PicardTools):

    # The name of the tool
    _tool_name = "MarkDuplicates"

    # The jar location
    _jar = "/opt/picard-tools-1.113/MarkDuplicates.jar"

    # The options
    _command = ("INPUT={input} O={output} METRICS_FILE={metrics}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":   GenericTool.INPUT,
                         "output":  GenericTool.OUTPUT,
                         "metrics": GenericTool.OUTPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "dedup"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
       """Initialize a MarkDuplicates instance."""
       pass

    def execute(self, options, out_dir=None, locally=True):
        """Mark duplicates."""
        # The metrics file
        if "output" not in options:
            m = "{}: no output file".format(self.__class__.__name__)
            raise ProgramError(m)
        metrics_file = re.sub(r"\.[sb]am$", ".dedup", options["output"])
        options["metrics"] = metrics_file
        super(MarkDuplicates, self).execute(options, out_dir, locally)
