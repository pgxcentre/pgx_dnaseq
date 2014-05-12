__all__ = ["RealignerTargetCreator", "IndelRealigner", "PrintReads",
           "BaseRecalibrator", "UnifiedGenotyper", "UnifiedGenotyper_Multi",
           "HaplotypeCaller", "HaplotypeCaller_Multi"]

import os
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
                "-o {output} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL}

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
                "-targetIntervals {interval_file} -o {output} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":         GenericTool.INPUT,
                         "output":        GenericTool.OUTPUT,
                         "reference":     GenericTool.INPUT,
                         "interval_file": GenericTool.INPUT,
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
                "-o {output} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "groups":    GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL}

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
                "-knownSites {dbsnp} -o {groups} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{groups}.out"
    _stderr = "{groups}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "groups":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "dbsnp":     GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL}

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


class UnifiedGenotyper(GATK):

    # The name of the tool
    _tool_name = "UnifiedGenotyper"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T UnifiedGenotyper -R {reference} -I {input} --dbsnp {dbsnp} "
                "-o {output} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "dbsnp":     GenericTool.INPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "unified_genotyper"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.vcf".format(_suffix), )

    def __init__(self):
        """Initialize a UnifiedGenotyper instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes a BAM and calls."""
        # First we index the input file
        if "input" not in options:
            m = "{}: no input file".format(self.__class__.__name__)
            raise ProgramError(m)
        IndexBam().execute({"input": options["input"]}, out_dir)

        # Then we create the MPILEUP file
        super(UnifiedGenotyper, self).execute(options, out_dir)


class UnifiedGenotyper_Multi(GATK):

    # The name of the tool
    _tool_name = "UnifiedGenotyper_Multi"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T UnifiedGenotyper -R {reference} --input_file {input} "
                "--dbsnp {dbsnp} -o {output} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "dbsnp":     GenericTool.INPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "unified_genotyper"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.vcf".format(_suffix), )

    # This tool needs multiple input
    _merge_all_inputs = True

    def __init__(self):
        """Initialize a UnifiedGenotyper_Multi instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes all the BAM files and calls."""
        # First we index all the input files
        if "inputs" not in options:
            m = "{}: no input files".format(self.__class__.__name__)
            raise ProgramError(m)
        for filename in options["inputs"]:
            IndexBam().execute({"input": filename}, out_dir)

        # We need to create the list of bam files
        list_filename = os.path.join(out_dir, "input_files.list")
        with open(list_filename, "w") as o_file:
            print("\n".join(options["inputs"]), file=o_file)

        # The input file is now the file containing the list
        options["input"] = list_filename

        # Then we call
        super(UnifiedGenotyper_Multi, self).execute(options, out_dir)


class HaplotypeCaller(GATK):

    # The name of the tool
    _tool_name = "HaplotypeCaller"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T HaplotypeCaller -R {reference} -I {input} --dbsnp {dbsnp} "
                "-o {output} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "dbsnp":     GenericTool.INPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "haplotype_caller"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.vcf".format(_suffix), )

    def __init__(self):
        """Initialize a HaplotypeCaller instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes a BAM and calls."""
        # First we index the input file
        if "input" not in options:
            m = "{}: no input file".format(self.__class__.__name__)
            raise ProgramError(m)
        IndexBam().execute({"input": options["input"]}, out_dir)

        # Then we create the MPILEUP file
        super(HaplotypeCaller, self).execute(options, out_dir)


class HaplotypeCaller_Multi(GATK):

    # The name of the tool
    _tool_name = "HaplotypeCaller_Multi"

    # The jar file
    _jar = "GenomeAnalysisTK.jar"

    # The options
    _command = ("-T HaplotypeCaller -R {reference} --input_file {input} "
                "--dbsnp {dbsnp} -o {output} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL,
                         "dbsnp":     GenericTool.INPUT}

    # The suffix that will be added just before the extension of the output file
    _suffix = "haplotype_caller"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.vcf".format(_suffix), )

    # This tool needs multiple input
    _merge_all_inputs = True

    def __init__(self):
        """Initialize a HaplotypeCaller instance."""
        pass

    def execute(self, options, out_dir=None):
        """Indexes all the BAM files and calls."""
        # First we index all the input files
        if "inputs" not in options:
            m = "{}: no input files".format(self.__class__.__name__)
            raise ProgramError(m)
        for filename in options["inputs"]:
            IndexBam().execute({"input": filename}, out_dir)

        # We need to create the list of bam files
        list_filename = os.path.join(out_dir, "input_files.list")
        with open(list_filename, "w") as o_file:
            print("\n".join(options["inputs"]), file=o_file)

        # The input file is now the file containing the list
        options["input"] = list_filename

        # Then we call
        super(HaplotypeCaller_Multi, self).execute(options, out_dir)
