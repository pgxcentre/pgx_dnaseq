#!/usr/bin/env python3


__author__ = "Louis-Philippe Lemieux Perreault and Abdellatif Daghrach"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault", "Abdellatif Daghrach",
               "Michal Blazejczyk"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


import os
import re
import sys
import logging
import argparse
from statistics import mean, stdev
from collections import defaultdict
from subprocess import Popen, PIPE, TimeoutExpired

import jinja2
import matplotlib.pyplot as plt
from pkg_resources import resource_filename

from pgx_dna_seq import ProgramError
from pgx_dna_seq.read_config import get_pipeline_steps


def main():
    """The main function."""
    # Creating the option parser
    desc = "Generate automatic report (version {}).".format(__version__)
    parser = argparse.ArgumentParser(description=desc)

    # Running the script
    try:
        # Getting and checking the options
        args = parse_args(parser)
        check_args(args)

        # Logging
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=[logging.StreamHandler()],
        )

        # Logging the options used
        log_options(args)

        # Gathering the samples
        sample_list = gather_samples(args.input)

        # Getting the pipeline steps
        steps = get_pipeline_steps(args.pipeline_config)

        # Printing the report
        print_report(sample_list, steps, args)

    except KeyboardInterrupt:
        logging.info("Cancelled by user")
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    except ProgramError as e:
        logging.error(e)
        parser.error(e.message)

    except NotImplementedError as e:
        logging.error(e)

    except Exception as e:
        logging.error(e)
        raise


def log_options(options):
    """Logs the options used."""
    logging.info("The following options were used:")
    for opt, value in vars(options).items():
        logging.info("  --{} {}".format(opt.replace("_", "-"), value))


def gather_samples(filename):
    """Gather the sample list."""
    all_samples = set()
    with open(filename, "r") as i_file:
        for line in i_file:
            samples = re.split(r"\s+", line.rstrip("\n"))
            samples = [os.path.basename(i) for i in samples]

            for sample in samples:
                if sample.endswith(".fastq") or sample.endswith(".fastq.gz"):
                    all_samples.add(re.search(r"(^\w+)_R[12]", sample).group(1))

                else:
                    all_samples.add(re.search(r"\w+", sample).group())

    logging.info("Found {:,d} samples".format(len(all_samples)))
    return sorted(all_samples)


def print_report(sample_list, pipeline_steps, options):
    """Creates the report."""
    # Creating the jinja2 environment
    jinja2_env = jinja2.Environment(
        block_start_string = '\BLOCK{',
        block_end_string = '}',
        variable_start_string = '\VAR{',
        variable_end_string = '}',
        comment_start_string = '\#{',
        comment_end_string = '}',
        line_statement_prefix = '%-',
        line_comment_prefix = '%#',
        trim_blocks = True,
        autoescape = False,
        loader=jinja2.PackageLoader("pgx_dna_seq", "report_templates")
    )

    # Getting the flowchart (either PNG or PDF only)
    flowchart = None
    if os.path.isfile("flowchart.pdf"):
        flowchart = "flowchart.pdf"
    elif os.path.isfile("flowchart.png"):
        flowchart = "flowchart.png"

    # Getting the report content
    report_content, metrics_from = construct_report_content(sample_list,
                                                            pipeline_steps,
                                                            jinja2_env)

    # Creating the data to put into the final report
    resource_prefix = "report_templates/images/"
    report_data = {
        "logo_small":     resource_filename("pgx_dna_seq", resource_prefix +
                                            "pgx_logo_small.png"),
        "logo":           resource_filename("pgx_dna_seq",
                                            resource_prefix + "pgx_logo.pdf"),
        "run_name":       sanitize_tex(options.run_name),
        "flowchart":      flowchart,
        "report_content": report_content,
        "summary":        sanitize_tex(get_pipeline_summary(pipeline_steps)),
        "metrics_from":   sanitize_tex(pretty_join(
            [r"\texttt{" + name + "}" for name in sorted(metrics_from)]
        )),
    }

    # The name of the TEX file
    report_filename = re.sub(r"\.pdf$", ".tex", options.output)

    # Creating the TEX using the main template
    logging.info("Compiling the report")
    try:
        template = jinja2_env.get_template("main_template.tex")
        with open(report_filename, "w") as o_file:
            print(template.render(**report_data), file=o_file)

    except FileNotFoundError:
        raise ProgramError("{}: cannot write file".format(report_filename))

    # The output directory
    out_dir = os.path.dirname(options.output)
    out_dir = out_dir if out_dir else "."

    # Compiling the LaTeX report (2 times)
    command = ["pdflatex", "-halt-on-error", "-output-directory", out_dir,
               report_filename]
    try:
        for i in range(2):
            # Executing the command
            proc = Popen(command, stdout=PIPE, stderr=PIPE)

            # Waiting for the process to terminate
            try:
                outs, errs = proc.communicate(timeout=60)

            except TimeoutExpired:
                # Killing the process
                proc.kill()
                m = ("something went wrong while compiling the report: try "
                     "running the following "
                     "command\n\t{}".format(" ".join(command)))
                raise ProgramError(m)

            # Getting the return code
            rc = proc.returncode
            if rc != 0:
                # There was a problem...
                m = ("problem while compiling the report: check log "
                     "'{}'".format(re.sub(r"\.pdf", ".log", options.output)))
                raise ProgramError(m)

    except FileNotFoundError:
        m = "'{}' is not installed".format(command[0])
        raise ProgramError(m)

    # Now, we want to delete the following extensions
    ext_to_delete = [".aux", ".lot", ".log", ".toc", ".out", ".lof", ".tex"]
    for ext in ext_to_delete:
        file_to_delete = re.sub(r"\.pdf", ext, options.output)
        if os.path.isfile(file_to_delete):
            os.remove(file_to_delete)


def construct_report_content(samples, pipeline_steps, jinja2_env):
    """Constructs the report content."""
    report_content = ""

    available_steps = set()
    final_data = defaultdict(dict)
    for i, (step, step_option) in enumerate(pipeline_steps):
        prefix = os.path.join("output",
                              "{:02d}_{}".format(i+1, step.get_tool_name()))

        if step.get_tool_name() == "HsMetrics":
            available_steps.add("HsMetrics")
            logging.info("Collecting HsMetrics")
            for sample in samples:
                sample_data = step.read_report(os.path.join(prefix, sample))

                # The required values and type
                req_values = {
                    "total_reads":           int,
                    "off_bait_bases":        int,
                    "mean_bait_coverage":    float,
                    "mean_target_coverage":  float,
                    "fold_enrichment":       float,
                    "zero_cvg_targets_pct":  float,
                    "pct_target_bases_2x":   float,
                    "pct_target_bases_10x":  float,
                    "pct_target_bases_20x":  float,
                    "pct_target_bases_30x":  float,
                    "pct_target_bases_40x":  float,
                    "pct_target_bases_50x":  float,
                    "pct_target_bases_100x": float,
                }

                # Checking the required columns
                final_data = gather_values(sample_data, req_values, final_data,
                                           sample, prefix)

        if step.get_tool_name() == "MarkDuplicates":
            available_steps.add("MarkDuplicates")
            logging.info("Collecting MarkDuplicates")
            for sample in samples:
                sample_data = step.read_report(os.path.join(prefix, sample))

                # The required values and type
                req_values = {
                    "read_pair_duplicates":         int,
                    "percent_duplication":          float,
                    "read_pair_optical_duplicates": int,
                }

                # Checking the required columns
                final_data = gather_values(sample_data, req_values, final_data,
                                           sample, prefix)

        if step.get_tool_name() == "InsertSize":
            available_steps.add("InsertSize")
            logging.info("Collecting InsertSize")
            for sample in samples:
                sample_data = step.read_report(os.path.join(prefix, sample))

                # The required values and type
                req_values = {
                    "mean_insert_size":   float,
                    "median_insert_size": float,
                    "standard_deviation": float,
                    "hist_figname":       str,
                }

                # Checking the required columns
                final_data = gather_values(sample_data, req_values, final_data,
                                           sample, prefix)

        if step.get_tool_name() == "ClipTrim":
            available_steps.add("ClipTrim")
            logging.info("Collecting ClipTrim")
            for sample in samples:
                sample_data = step.read_report(os.path.join(prefix, sample))

                # The required values and type
                req_values = {
                    "total_reads_before_trim":   int,
                    "nb_short_reads_after_trim": int,
                    "nb_trimmed_r1":             int,
                    "nb_trimmed_r2":             int,
                }

                final_data = gather_values(sample_data, req_values, final_data,
                                           sample, prefix)

        if step.get_tool_name() == "CoverageGraph":
            available_steps.add("CoverageGraph")
            logging.info("Collecting CoverageGraph")
            for sample in samples:
                sample_data = step.read_report(os.path.join(prefix, sample))

                # The required values and type
                req_values = {
                    "coverage_figname": str,
                }

                final_data = gather_values(sample_data, req_values, final_data,
                                           sample, prefix)

        if step.get_tool_name() == "CoverageGraph_Multi":
            available_steps.add("CoverageGraph_Multi")
            logging.info("Collecting CoverageGraph_Multi")
            data = step.read_report(os.path.join(prefix, "all_samples"))

            # The required values and type
            req_values = {
                "coverage_multi_figname": str,
            }

            final_data = gather_values(data, req_values, final_data, sample,
                                       prefix)

        if step.get_tool_name() == "ReadQualityGraph":
            available_steps.add("ReadQualityGraph")
            logging.info("Collectiong ReadQualityGraph")
            for sample in samples:
                sample_data = step.read_report(os.path.join(prefix, sample))

                # The required values and type
                req_values = {
                    "read_qual_figname": str,
                }
                final_data = gather_values(sample_data, req_values, final_data,
                                           sample, prefix)

    sample_summary = defaultdict(list)
    final_tables = []
    sample_order = sorted(final_data.keys())
    for sample in sample_order:
        sample_data = final_data[sample]
        sample_tables = []

        # Clipping and trimming
        if "ClipTrim" in available_steps:
            total_reads = sample_data["total_reads_before_trim"]
            too_short = sample_data["nb_short_reads_after_trim"]
            trimmed_r1 = sample_data["nb_trimmed_r1"]
            trimmed_r2 = sample_data["nb_trimmed_r2"]

            # The table
            sample_tables.append({
                "name":   "Clipping/Trimming",
                "format": "lr",
                "data":   [["Total reads", "{:,d}".format(total_reads)],
                           ["Too short after clip", "{:,d}".format(too_short)],
                           ["Trimmed R1", "{:,d}".format(trimmed_r1)],
                           ["Trimmed R2", "{:,d}".format(trimmed_r2)]],
            })

            # Saving for the sample_summary
            sample_summary["total_reads_before_trim"].append(total_reads)
            sample_summary["nb_short_reads_after_trim"].append(too_short)
            sample_summary["nb_trimmed_r1"].append(trimmed_r1)
            sample_summary["nb_trimmed_r2"].append(trimmed_r2)

        # Mark duplicates
        if "MarkDuplicates" in available_steps:
            rp_dup = sample_data["read_pair_duplicates"]
            pc_dup = sample_data["percent_duplication"] * 100
            optical_dup = sample_data["read_pair_optical_duplicates"]
            optical_dup = optical_dup / rp_dup * 100
            pcr_dup = 100 - optical_dup

            # The table
            sample_tables.append({
                "name":   "Duplicated Reads",
                "format": "lr",
                "data":   [
                    ["Total duplicates", "{:,d}".format(rp_dup)],
                    ["Duplicated percentage", r"{:,.2f}\%".format(pc_dup)],
                    ["Optical duplicate", r"{:,.2f}\%".format(optical_dup)],
                    ["PCR duplicate", r"{:,.2f}\%".format(pcr_dup)]],
            })

            # Saving for the sample summary
            sample_summary["read_pair_duplicates"].append(rp_dup)
            sample_summary["percent_duplication"].append(pc_dup)
            sample_summary["read_pair_optical_duplicates"].append(optical_dup)
            sample_summary["read_pair_pcr_duplicates"].append(pcr_dup)

        # HsMetrics
        if "HsMetrics" in available_steps:
            total_reads = sample_data["total_reads"]
            off_bait = sample_data["off_bait_bases"]
            mean_bait = sample_data["mean_bait_coverage"]
            mean_target = sample_data["mean_target_coverage"]
            fold = sample_data["fold_enrichment"]
            cvg_targets = sample_data["zero_cvg_targets_pct"]
            pct_target_2x = sample_data["pct_target_bases_2x"]
            pct_target_10x = sample_data["pct_target_bases_10x"]
            pct_target_20x = sample_data["pct_target_bases_20x"]
            pct_target_30x = sample_data["pct_target_bases_30x"]
            pct_target_40x = sample_data["pct_target_bases_40x"]
            pct_target_50x = sample_data["pct_target_bases_50x"]
            pct_target_100x = sample_data["pct_target_bases_100x"]

            # The table
            sample_tables.append({
                "name":   "HS Metrics",
                "format": "lr",
                "data":   [
                    ["Total reads", "{:,d}".format(total_reads)],
                    ["Off bait", "{:,d}".format(off_bait)],
                    ["Mean bait coverage", "{:,.2f}".format(mean_bait)],
                    ["Mean target coverage", "{:,.2f}".format(mean_target)],
                    ["Fold enrichment", "{:,.2f}".format(fold)],
                    ["Zero CVG target", r"{:,.2f}\%".format(cvg_targets)],
                    ["Target bases (2x)", r"{:,.2f}\%".format(pct_target_2x)],
                    ["Target bases (10x)", r"{:,.2f}\%".format(pct_target_10x)],
                    ["Target bases (20x)", r"{:,.2f}\%".format(pct_target_20x)],
                    ["Target bases (30x)", r"{:,.2f}\%".format(pct_target_30x)],
                    ["Target bases (40x)", r"{:,.2f}\%".format(pct_target_40x)],
                    ["Target bases (50x)", r"{:,.2f}\%".format(pct_target_50x)],
                    ["Target bases (100x)", r"{:,.2f}\%".format(pct_target_100x)]],
            })

            # Saving for the sample summary
            sample_summary["total_reads"].append(total_reads)
            sample_summary["off_bait_bases"].append(off_bait)
            sample_summary["mean_bait_coverage"].append(mean_bait)
            sample_summary["mean_target_coverage"].append(mean_target)
            sample_summary["fold_enrichment"].append(fold)
            sample_summary["zero_cvg_targets_pct"].append(cvg_targets)
            sample_summary["pct_target_bases_2x"].append(pct_target_2x)
            sample_summary["pct_target_bases_10x"].append(pct_target_10x)
            sample_summary["pct_target_bases_20x"].append(pct_target_20x)
            sample_summary["pct_target_bases_30x"].append(pct_target_30x)
            sample_summary["pct_target_bases_40x"].append(pct_target_40x)
            sample_summary["pct_target_bases_50x"].append(pct_target_50x)
            sample_summary["pct_target_bases_100x"].append(pct_target_100x)

        # InsertSize
        size_hist = None
        if "InsertSize" in available_steps:
            mean_size = sample_data["mean_insert_size"]
            median_size = sample_data["median_insert_size"]
            std_size = sample_data["standard_deviation"]

            # The table
            sample_tables.append({
                "name":   "Insert Size",
                "format": "lr",
                "data":   [["Mean", "{:,.2f}".format(mean_size)],
                           ["Median", "{:,.2f}".format(median_size)],
                           ["Standard deviation", "{:,.2f}".format(std_size)]],
            })

            # The histogram
            size_hist = sample_data["hist_figname"]

            # Saving for the sample summary
            sample_summary["mean_insert_size"].append(mean_size)
            sample_summary["median_insert_size"].append(median_size)
            sample_summary["standard_deviation"].append(std_size)

        # The coverage plot (per sample)
        cov_plot = None
        if "CoverageGraph" in available_steps:
            cov_plot = sample_data["coverage_figname"]

        # The read quality graph (per sample)
        qual_plot = None
        if "ReadQualityGraph" in available_steps:
            qual_plot = sample_data["read_qual_figname"]

        # Saving tables
        final_tables.append({
            "name": sample,
            "tables": sample_tables,
            "plot_1": size_hist,
            "plot_2": qual_plot,
            "plot_3": cov_plot
        })

    # The summary plots
    cov_multi_plot = None
    if "CoverageGraph_Multi" in available_steps:
        cov_multi_plot = sample_data["coverage_multi_figname"]

    # Getting the report content
    template = jinja2_env.get_template("data_template.tex")
    report_content = generate_sample_summary(sample_summary, sample_order,
                                             available_steps, cov_multi_plot,
                                             template)
    for data in final_tables:
        # Acquiring the data
        section_name = data["name"]
        section_tables = data["tables"]
        plot_1 = data["plot_1"]
        plot_2 = data["plot_2"]
        plot_3 = data["plot_3"]

        # If there is only one None plot, it should be the second one...
        first_plot = None
        second_plot = None
        for plot in [plot_1, plot_2, plot_3]:
            if plot is not None:
                if first_plot is None:
                    first_plot = plot
                elif second_plot is None:
                    second_plot = plot

        # The template data
        report_data = {
            "section_name": section_name,
            "tables":       section_tables,
            "first_plot":   first_plot,
            "second_plot":  second_plot,
        }

        # Generating the template
        report_content += template.render(**report_data)

    return report_content, available_steps


def generate_sample_summary(sample_values, sample_order, steps, cov_multi,
                            template):
    """Generates the sample summary."""
    if len(sample_order) < 2:
        return ""

    summary_tables = []
    # Clipping and trimming
    if "ClipTrim" in steps:
        total_reads = sample_values["total_reads_before_trim"]
        too_short = sample_values["nb_short_reads_after_trim"]
        trimmed_r1 = sample_values["nb_trimmed_r1"]
        trimmed_r2 = sample_values["nb_trimmed_r2"]

        # The table
        summary_tables.append({
            "name":   "Clipping/Trimming",
            "format": "lrl",
            "data":   [[
                "Total reads",
                "{:,.2f}".format(mean(total_reads)),
                r"$\pm$ {:,.2f}".format(stdev(total_reads)),
            ],
            [
                "Too short after clip",
                "{:,.2f}".format(mean(too_short)),
                r"$\pm$ {:,.2f}".format(stdev(too_short)),
            ],
            [
                "Trimmed R1",
                "{:,.2f}".format(mean(trimmed_r1)),
                r"$\pm$ {:,.2f}".format(stdev(trimmed_r1)),
            ],
            [
                "Trimmed R2",
                "{:,.2f}".format(mean(trimmed_r2)),
                r"$\pm$ {:,.2f}".format(stdev(trimmed_r2)),
            ]],
        })

    # Mark duplicates
    if "MarkDuplicates" in steps:
        rp_dup = sample_values["read_pair_duplicates"]
        pc_dup = sample_values["percent_duplication"]
        optical_dup = sample_values["read_pair_optical_duplicates"]
        pcr_dup = sample_values["read_pair_pcr_duplicates"]

        # The table
        summary_tables.append({
            "name":   "Duplicated Reads",
            "format": "lrl",
            "data":   [[
                "Total duplicates",
                "{:,.2f}".format(mean(rp_dup)),
                r"$\pm$ {:,.2f}".format(stdev(rp_dup)),
            ],
            [
                "Duplicated percentage",
                r"{:,.2f}\%".format(mean(pc_dup)),
                r"$\pm$ {:,.2f}".format(stdev(pc_dup)),
            ],
            [
                "Optical duplicate",
                r"{:,.2f}\%".format(mean(optical_dup)),
                r"$\pm$ {:,.2f}".format(stdev(optical_dup)),
            ],
            [
                "PCR duplicate",
                r"{:,.2f}\%".format(mean(pcr_dup)),
                r"$\pm$ {:,.2f}".format(stdev(pcr_dup)),
            ]],
        })

    # HsMetrics
    if "HsMetrics" in steps:
        total_reads = sample_values["total_reads"]
        off_bait = sample_values["off_bait_bases"]
        mean_bait = sample_values["mean_bait_coverage"]
        mean_target = sample_values["mean_target_coverage"]
        fold = sample_values["fold_enrichment"]
        cvg_targets = sample_values["zero_cvg_targets_pct"]
        pct_target_2x = sample_values["pct_target_bases_2x"]
        pct_target_10x = sample_values["pct_target_bases_10x"]
        pct_target_20x = sample_values["pct_target_bases_20x"]
        pct_target_30x = sample_values["pct_target_bases_30x"]
        pct_target_40x = sample_values["pct_target_bases_40x"]
        pct_target_50x = sample_values["pct_target_bases_50x"]
        pct_target_100x = sample_values["pct_target_bases_100x"]

        # The table
        summary_tables.append({
            "name":   "HS Metrics",
            "format": "lrl",
            "data":   [[
                "Total reads",
                "{:,.2f}".format(mean(total_reads)),
                r"$\pm$ {:,.2f}".format(stdev(total_reads)),
            ],
            [
                "Off bait",
                "{:,.2f}".format(mean(off_bait)),
                r"$\pm$ {:,.2f}".format(stdev(off_bait)),
            ],
            [
                "Mean bait coverage",
                "{:,.2f}".format(mean(mean_bait)),
                r"$\pm$ {:,.2f}".format(stdev(mean_bait)),
            ],
            [
                "Mean target coverage",
                "{:,.2f}".format(mean(mean_target)),
                r"$\pm$ {:,.2f}".format(stdev(mean_target)),
            ],
            [
                "Fold enrichment",
                "{:,.2f}".format(mean(fold)),
                r"$\pm$ {:,.2f}".format(stdev(fold)),
            ],
            [
                "Zero CVG target",
                r"{:,.2f}\%".format(mean(cvg_targets)),
                r"$\pm$ {:,.2f}".format(stdev(cvg_targets)),
            ],
            [
                "Target bases (2x)",
                r"{:,.2f}\%".format(mean(pct_target_2x)),
                r"$\pm$ {:,.2f}".format(stdev(pct_target_2x)),
            ],
            [
                "Target bases (10x)",
                r"{:,.2f}\%".format(mean(pct_target_10x)),
                r"$\pm$ {:,.2f}".format(stdev(pct_target_10x)),
            ],
            [
                "Target bases (20x)",
                r"{:,.2f}\%".format(mean(pct_target_20x)),
                r"$\pm$ {:,.2f}".format(stdev(pct_target_20x)),
            ],
            [
                "Target bases (30x)",
                r"{:,.2f}\%".format(mean(pct_target_30x)),
                r"$\pm$ {:,.2f}".format(stdev(pct_target_30x)),
            ],
            [
                "Target bases (40x)",
                r"{:,.2f}\%".format(mean(pct_target_40x)),
                r"$\pm$ {:,.2f}".format(stdev(pct_target_40x)),
            ],
            [
                "Target bases (50x)",
                r"{:,.2f}\%".format(mean(pct_target_50x)),
                r"$\pm$ {:,.2f}".format(stdev(pct_target_50x)),
            ],
            [
                "Target bases (100x)",
                r"{:,.2f}\%".format(mean(pct_target_100x)),
                r"$\pm$ {:,.2f}".format(stdev(pct_target_100x)),
            ]],
        })

    # InsertSize
    size_hist = None
    if "InsertSize" in steps:
        mean_size = sample_values["mean_insert_size"]
        median_size = sample_values["median_insert_size"]
        std_size = sample_values["standard_deviation"]

        # The figure and axe
        figure, axe = plt.subplots(1, 1, figsize=(12, 6))
        axe.bar(range(len(mean_size)), mean_size,
                align="center", color="#0099CC", edgecolor="#0099CC", lw=0)

        # The labels
        axe.set_title("All samples mean insert size", weight="bold",
                      fontsize=12)
        axe.set_ylabel("Mean Size", weight="bold", fontsize=10)

        # The x tick labels
        axe.set_xticks(range(len(mean_size)))
        axe.set_xticklabels(sample_order)

        # The tick labels fontsize
        axe.tick_params(axis='both', which='major', labelsize=9)

        # Saving the figure
        size_hist = "{}.pdf".format(os.path.join("output", "mean_insert_size"))
        plt.savefig(size_hist, figure=figure, bbox_inches="tight")
        plt.close(figure)

        # The table
        summary_tables.append({
            "name":   "Insert Size",
            "format": "lrl",
            "data":   [[
                "Mean",
                "{:,.2f}".format(mean(mean_size)),
                r"$\pm$ {:,.2f}".format(stdev(mean_size)),
            ],
            [
                "Median",
                "{:,.2f}".format(mean(median_size)),
                r"$\pm$ {:,.2f}".format(stdev(median_size)),
            ],
            [
                "Standard deviation",
                "{:,.2f}".format(mean(std_size)),
                r"$\pm$ {:,.2f}".format(stdev(std_size)),
            ]],
        })

    # The plots
    first_plot = size_hist
    second_plot = cov_multi
    if first_plot is None and second_plot is not None:
        first_plot, second_plot = second_plot, first_plot

    # The template data
    report_data = {
        "section_name": "All samples",
        "tables":       summary_tables,
        "first_plot":   first_plot,
        "second_plot":  second_plot,
    }

    # Generating the template
    return template.render(**report_data)


def gather_values(values, required_values, final_values, sample, prefix):
    """Gather all the required values."""
    for name in required_values.keys():
        if name in final_values[sample]:
            logging.warning("{}: duplicated values (only the most recent will "
                            "be kept".format(name))
        if name not in values:
            raise ProgramError("ClipTrim: {}: no value named "
                                "{}".format(prefix, name))
        final_values[sample][name] = required_values[name](values[name])

    return final_values


def get_pipeline_summary(steps):
    """Gets the pipeline summary."""
    steps = [r"\texttt{" + i.get_tool_name() + "}" for i, j in steps]
    return pretty_join(steps)


def pretty_join(items):
    """Pretty join a list for English."""
    return ", ".join(items[:-1]) + " and " + items[-1]


def sanitize_tex(text):
    """Sanitize TeX text."""
    for char in ["_", "#", "%"]:
        text = text.replace(char, r"\{}".format(char))

    return text


def check_args(args):
    """Checks the arguments and options."""
    # Checking that the configuration file exists
    if not os.path.isfile(args.pipeline_config):
        m = "{}: no such file".format(args.pipeline_config)
        raise ProgramError(m)

    # Checking that the file containing all input files exists
    if not os.path.isfile(args.input):
        m = "{}: no such file".format(args.input)
        raise ProgramError(m)

    # Checking the output report
    if not args.output.endswith(".pdf"):
        args.output += ".pdf"
   
    return True


def parse_args(parser):
    """Parses the command line options and arguments."""
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s version {}".format(__version__))
    parser.add_argument("--debug", action="store_true",
                        help="Set the log level to debug.")

    group = parser.add_argument_group("General Options")
    group.add_argument("-i", "--input", type=str, metavar="FILE",
                       default="input_files.txt",
                       help=("A file containing the pipeline input files (one "
                             "sample per line, one or more file per sample "
                             "[%(default)s]"))
    group.add_argument("-p", "--pipeline-config", type=str, metavar="FILE",
                       default="pipeline.conf",
                       help="The pipeline configuration file [%(default)s]")

    group = parser.add_argument_group("Report Options")
    group.add_argument("-r", "--run-name", type=str, metavar="STRING",
                       default="Sequencing_Run",
                       help="The Sequencing run name [%(default)s]")
    group.add_argument("-o", "--output", type=str, metavar="FILE",
                       default="dna_seq_report.pdf",
                       help="The name of the final report [%(default)s]")

    return parser.parse_args()


if __name__ == "__main__":
    main()

