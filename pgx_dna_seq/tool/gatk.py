__all__ = ["RealignerTargetCreator", "IndelRealigner", "PrintReads",
           "BaseRecalibrator"]

import re

from pgx_dna_seq import ProgramError
from pgx_dna_seq.tool.java import JAR
from pgx_dna_seq.tool import GenericTool
from pgx_dna_seq.tool.samtools import IndexBam


class GATK(JAR):

    # The version of the tool
    _version = "3.1-1"

    # The jar file default location
    _jar_location = "/opt/GenomeAnalysisTK-3.1-1"

    def __init__(self):
        """Initialize a PicardTools instance."""
        pass


class RealignerTargetCreator(GATK):

    # The name of the tool
    _tool_name = "RealignerTargetCreator"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T RealignerTargetCreator -I {input} -R {reference} "
                "-o {output} -dt {dt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "dt":        GenericTool.REQUIREMENT}

    # The suffix that will be added just before the extension of the output file
    _suffix = None

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".intervals", )

    def __init__(self):
       """Initialize a _Realign instance."""
       pass


class IndelRealigner(GATK):

    # The name of the tool
    _tool_name = "IndelRealigner"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T IndelRealigner -R {reference} -I {input} "
                "-targetIntervals {interval_file} -o {output} -dt {dt} "
                "{other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":         GenericTool.INPUT,
                         "output":        GenericTool.OUTPUT,
                         "reference":     GenericTool.INPUT,
                         "interval_file": GenericTool.INPUT,
                         "dt":            GenericTool.REQUIREMENT,
                         "other_opt":     GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
    _suffix = "realign_gatk"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
        """Initialize a IndelRealigner instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes a BAM and realign it."""
        # First we index the input file
        if "input" not in options:
            m = "{}: no input file".format(self.__class__.__name__)
            raise ProgramError(m)
        IndexBam().execute({"input": options["input"]}, out_dir)

        # Then we create the intervals for realigning
        intervals_opt = {}
        try:
            # The options
            intervals_opt["input"] = options["input"]
            intervals_opt["reference"] = options["reference"]
            intervals_opt["dt"] = options["dt"]
            if "java_memory" in options:
                intervals_opt["java_memory"] = options["java_memory"]

            # The intervals filename
            intervals_filename = re.sub(r"\.[sb]am$", ".intervals",
                                        options["output"])
            intervals_opt["output"] = intervals_filename
            options["interval_file"] = intervals_filename
        except KeyError as e:
            m = "{}: missing option {}".format(self.__class__.__name__, e)
            raise ProgramError(m)

        # Executing the intervals
        RealignerTargetCreator().execute(intervals_opt, out_dir)

        # Executing the realignment
        super(IndelRealigner, self).execute(options, out_dir)


class PrintReads(GATK):

    # The name of the tool
    _tool_name = "PrintReads"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T PrintReads -I {input} -R {reference} -BQSR {groups} "
                "-o {output} -dt {dt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "groups":    GenericTool.INPUT,
                         "dt":        GenericTool.REQUIREMENT}

    # The suffix that will be added just before the extension of the output file
    _suffix = None

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".grp", )

    def __init__(self):
       """Initialize a PrintReads instance."""
       pass


class BaseRecalibrator(GATK):

    # The name of the tool
    _tool_name = "BaseRecalibrator"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T BaseRecalibrator -R {reference} -I {input} "
                "-knownSites {dbSNP_known_sites} -o {groups} -dt {dt}")

    # The STDOUT and STDERR
    _stdout = "{groups}.out"
    _stderr = "{groups}.err"

    # The description of the required options
    _required_options = {"input":             GenericTool.INPUT,
                         "groups":            GenericTool.OUTPUT,
                         "reference":         GenericTool.INPUT,
                         "dbSNP_known_sites": GenericTool.INPUT,
                         "dt":                GenericTool.REQUIREMENT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "base_recal"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
        """Initialize a BaseRecalibrator instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes a BAM and recalibrate it."""
        # First we index the input file
        if "input" not in options:
            m = "{}: no input file".format(self.__class__.__name__)
            raise ProgramError(m)
        IndexBam().execute({"input": options["input"]}, out_dir)

        # Then we create the groups for recalibration
        if "output" not in options:
            m = "{}: no output file".format(self.__class__.__name__)
            raise ProgramError(m)
        groups_file = re.sub(r"\.[sb]am$", ".grp", options["output"])
        options["groups"] = groups_file
        super(BaseRecalibrator, self).execute(options, out_dir)

        # Printing the reads
        PrintReads().execute(options, out_dir)
