
# This file is part of pgx_dnaseq
#
# This work is licensed under The MIT License (MIT). To view a copy of this
# license, visit http://opensource.org/licenses/MIT


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "MIT"


import os
import re
from glob import glob

from .java import JAR
from . import GenericTool


__all__ = ["SortSam", "AddRG", "MarkDuplicates", "HsMetrics" , "InsertSize"]


class PicardTools(JAR):

    # The version of the tool
    _version = "1.127"

    # The jar file default location
    _jar_location = "/opt/picard-tools-1.127"
    _jar = "picard.jar"

    def __init__(self):
        """Initialize a PicardTools instance."""
        pass


class SortSam(PicardTools):

    # The name of the tool
    _tool_name = "SortSam"

    # The options
    _command = ("SortSam INPUT={input} OUTPUT={output} SORT_ORDER={sort_order} "
                "{other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":      GenericTool.INPUT,
                         "output":     GenericTool.OUTPUT,
                         "sort_order": GenericTool.REQUIREMENT,
                         "other_opt":  GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
    _suffix = "sorted"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
       """Initialize a SortSam instance."""
       pass


class HsMetrics(PicardTools):

    # The name of the tool
    _tool_name = "HsMetrics"

    # The options
    _command = ("CalculateHsMetrics INPUT={input} OUTPUT={output} "
                "REFERENCE_SEQUENCE={reference} BAIT_INTERVALS={baits} "
                "TARGET_INTERVALS={targets} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "baits":     GenericTool.INPUT,
                         "targets":   GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
    _suffix = "hsmetrics"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}".format(_suffix), )

    # This tool does not produce usable data...
    _produce_data = False
    
    def __init__(self):
       """Initialize a HSmetrics instance."""
       pass

    def read_report(self, prefix):
        """Reads a HsMetrics report file."""
        # Getting the report file name
        filename = glob("{}*.{}".format(prefix, self._suffix))
        assert len(filename) == 1
        filename = filename[0]

        result = {}
        with open(filename, "r") as i_file:
            # Reading until the header
            header = i_file.readline()
            while not header.startswith("## METRICS"):
                header = i_file.readline()

            # The two rows (header and data)
            header = i_file.readline().rstrip("\n").split("\t")
            data = i_file.readline().rstrip("\n").split("\t")

            for name, value in zip(header, data):
                result[name.lower()] = value

        return result


class InsertSize(PicardTools):
    # The name of the tool
    _tool_name = "InsertSize"

    # The options
    _command = ("CollectInsertSizeMetrics INPUT={input} OUTPUT={output} "
                "REFERENCE_SEQUENCE={reference} HISTOGRAM_FILE={hist_file} "
                "{other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "hist_file": GenericTool.OUTPUT,
                         "reference": GenericTool.INPUT,
                         "other_opt": GenericTool.OPTIONAL}
    
    # The suffix that will be added just before the extension of the output file
    _suffix = "insertsize"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}".format(_suffix), )
    
    # This tool does not produce usable data...
    _produce_data = False
    
    def __init__(self):
        """Initialize a HSmetrics instance."""
    pass
    
    def execute(self, options, out_dir=None):
        """InsertSize."""
        # The plot file
        if "output" not in options:
            m = "{}: no output file".format(self.__class__.__name__)
            raise ProgramError(m)
        plot_file =re.sub(r"\..*",".png",options["output"])
        options["hist_file"]=plot_file
        super(InsertSize,self).execute(options,out_dir)

    def read_report(self, prefix):
        """Reads a InsertSize report file."""
        import pandas as pd
        import matplotlib.pyplot as plt

        # Getting the report file name
        filename = glob("{}*.{}".format(prefix, self._suffix))
        assert len(filename) == 1
        filename = filename[0]

        result = {}
        with open(filename, "r") as i_file:
            # Reading until the header
            header = i_file.readline()
            while header.startswith("#") or header == "\n":
                header = i_file.readline()

            # The two rows (header and data)
            header = header.rstrip("\n").split("\t")
            data = i_file.readline().rstrip("\n").split("\t")

            for name, value in zip(header, data):
                result[name.lower()] = value

            # Then, reading until the histogram
            line = i_file.readline()
            while not line.startswith("## HISTOGRAM"):
                line = i_file.readline()

            # Plotting the barplot
            hist_data = pd.read_csv(i_file, sep="\t")

            # The figure and axe
            figure, axe = plt.subplots(1, 1, figsize=(12, 5))
            axe.bar(hist_data.insert_size, hist_data["All_Reads.fr_count"],
                    align="center", color="#0099CC", edgecolor="#0099CC", lw=0)

            # The labels
            axe.set_title("{} - Insert Size".format(os.path.basename(prefix)),
                          weight="bold", fontsize=12)
            axe.set_xlabel("Size", weight="bold", fontsize=10)
            axe.set_ylabel("Count", weight="bold", fontsize=10)

            # The tick labels fontsize
            axe.tick_params(axis='both', which='major', labelsize=8)

            # Saving the figure
            figname = "{}_{}.pdf".format(prefix, self._suffix)
            plt.savefig(figname, figure=figure, bbox_inches="tight")
            plt.close(figure)

            result["hist_figname"] = figname

        return result


class AddRG(PicardTools):

    # The name of the tool
    _tool_name = "AddRG"

    # The options
    _command = ("AddOrReplaceReadGroups I={input} O={output} RGDS={rgds} "
                "RGPL={rgpl} RGPU={rgpu} RGSM={sample_id} RGCN={rgcn} "
                "RGLB={rglb}")

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

    # The options
    _command = ("MarkDuplicates INPUT={input} O={output} "
                "METRICS_FILE={metrics} {other_opt}")

    # The STDOUT and STDERR
    _stdout = "{output}.out"
    _stderr = "{output}.err"

    # The description of the required options
    _required_options = {"input":     GenericTool.INPUT,
                         "output":    GenericTool.OUTPUT,
                         "metrics":   GenericTool.OUTPUT,
                         "other_opt": GenericTool.OPTIONAL}

    # The suffix that will be added just before the extension of the output file
    _suffix = "dedup"

    # The input and output type
    _input_type = (r"\.(\S+\.)?[sb]am$", )
    _output_type = (".{}.bam".format(_suffix), )

    def __init__(self):
       """Initialize a MarkDuplicates instance."""
       pass

    def execute(self, options, out_dir=None):
        """Mark duplicates."""
        # The metrics file
        if "output" not in options:
            m = "{}: no output file".format(self.__class__.__name__)
            raise ProgramError(m)
        metrics_file = re.sub(r"\.[sb]am$", ".dedup", options["output"])
        options["metrics"] = metrics_file
        super(MarkDuplicates, self).execute(options, out_dir)

    def read_report(self, prefix):
        """Reads a MarkDuplicates report file."""
        # Getting the report file name
        filename = glob("{}*.{}".format(prefix, self._suffix))
        assert len(filename) == 1
        filename = filename[0]

        result = {}
        with open(filename, "r") as i_file:
            # Reading until the header
            header = i_file.readline()
            while not header.startswith("## METRICS"):
                header = i_file.readline()

            # The two rows (header and data)
            header = i_file.readline().rstrip("\n").split("\t")
            data = i_file.readline().rstrip("\n").split("\t")

            for name, value in zip(header, data):
                result[name.lower()] = value

        return result

