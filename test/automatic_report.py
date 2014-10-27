#!/usr/bin/env python3

import os
import shutil
import re
import sys
import argparse
import glob

from datetime import date

import numpy as np


from ruffus import originate , formatter

import pgx_dna_seq

from pgx_dna_seq import ProgramError
from pgx_dna_seq.read_config import read_config_file, get_pipeline_steps

prog_version = "0.1"

# check input files
def check_input_files(filename):
    """Checks the input file names."""
    split_re = re.compile(r"\s+")
    input_filenames = None
    sample_list = list ()
    with open(filename, "r") as i_file:
        input_filenames = [re.split(r"\s+", i.rstrip("\r\n"))
        for i in i_file.readlines()]
    # Checking that all those files exists
    for sample_files in input_filenames:
        sample_list.append(sample_files[0].replace("data/","").replace("_R1.fastq.gz","")) 
        for input_filename in sample_files:
            if not os.path.isfile(input_filename):
                m = "{}: no such file".format(input_filename)
                raise ProgramError(m)
    
    return input_filenames,sample_list

# Read pipeline config file:

def check_args(args):
    """Checks the arguments and options.

    :param args: an object containing the options and arguments of the program.

    :type args: :py:class:`argparse.Namespace`

    ::returns: ``True`` if everything was OK.

    If there is a problem with an option, an exception is raised using the
    :py:class:`ProgramError` class, a message is printed to the
    :class:`sys.stderr` and the program exits with error code 1.

    """
    # Checking that the configuration file exists
    if not os.path.isfile(args.pipeline_config):
        m = "{}: no such file".format(args.pipeline_config)
        raise ProgramError(m)

    # Checking that the file containing all input files exists
    if not os.path.isfile(args.input):
        m = "{}: no such file".format(args.input)
        raise ProgramError(m)
   

    return True


def parse_args():
    """Parses the command line options and arguments.

    :returns: A :py:class:`argparse.Namespace` object created by the
              :py:mod:`argparse` module. It contains the values of the different
              options.

    ===============   =======  ================================================
        Options        Type                      Description
    ===============   =======  ================================================
    ===============   =======  ================================================

    .. note::
        No option check is done here (except for the one automatically done by
        :py:mod:`argparse`). Those need to be done elsewhere (see
        :py:func:`checkArgs`).

    """
    return parser.parse_args()


def print_report(what_to_run,sample_list):
    rmlist = ['report.toc','report.lof','report.out','report.log','report.aux']
    pages =list ()
    for sample in  (sample_list):   

        for job_index, (job, job_options) in enumerate(what_to_run):
            output_dir = os.path.join("output","{:02d}_{}".format(job_index +1,job.get_tool_name()))

        #Hsmetrics
            if job.get_tool_name() ==  'HsMetrics' :
                    hsmetrics_file=glob.glob(output_dir +'/'+sample+'*.hsmetrics')[0]
                    myfile=open (hsmetrics_file)
                    data= list ()
                    for i, line in enumerate (myfile):
                        if i>6:
                            data.append(line.replace("_","-").rstrip("\r\n").split("\t"))
                            Total_Reads=data[0][5]
                            Off_Bait=data[0][15]                
                            Mean_Bait_Coverage=data[0][20]      
                            Fold_Enrichement=data[0][24]        
                            Mean_Target_Coverage=data[0][21]    
                            Zero_CVG_Target=data[0][25]         
                            PCT_TARGET_BASES_2X=data[0][27]     
                            PCT_TARGET_BASES_10X=data[0][28]    
                            PCT_TARGET_BASES_20X=data[0][29]    
                            PCT_TARGET_BASES_30X=data[0][30]    
                            PCT_TARGET_BASES_40X=data[0][31]    
                            PCT_TARGET_BASES_50X=data[0][32]    
                            PCT_TARGET_BASES_100X=data[0][33]   

         #Duplicate 
            if job.get_tool_name() ==  'MarkDuplicates' :
                    dedup_file=glob.glob(output_dir +'/'+sample+'*.dedup')[0]
                    myfile=open (dedup_file)
                    data= list ()
                    for i, line in enumerate (myfile):
                        if i>6:
                            data.append(line.replace("_","-").rstrip("\r\n").split("\t"))
                            Total_Duplicate=data[0][5]
                            Optical_Duplicate=str (int(data[0][6]) * 100 / int(Total_Duplicate))
                            PCR_Duplicate=str(100-float(Optical_Duplicate))
                            Duplicate_Percentage=str(float(data[0][7]) * 100 )
         #Coverage_graph
            if job.get_tool_name() ==  'CoverageGraph' :
                    print (sample)
                    coverage_graph = glob.glob(output_dir +'/'+sample+'*.png')[0]
                    coverage_graph=coverage_graph.replace(".png",".png").replace("/"+sample,"/"+sample)
                    shutil.copy(coverage_graph, output_dir +'/'+sample+'.png')
                    coverage_graph=output_dir +'/'+sample+'.png'
                    rmlist.append(coverage_graph)
        #CoverageGraph_Multi
            if (job.get_tool_name() ==  'CoverageGraph_Multi' and 'coverage_graph_multi'  not  in locals()):
                    coverage_graph_multi = glob.glob(output_dir+'/all_samples*.png')[0]
                    shutil.copy(coverage_graph_multi, output_dir+'/all_samples.png')
                    coverage_graph_multi=output_dir +'/all_samples.png'
                    rmlist.append(coverage_graph_multi)
         #ClipTrim 
            if job.get_tool_name() ==  'ClipTrim' :
                    clip_file=output_dir +'/'+sample+'.out'
                    myfile=open (clip_file)
                    data= list ()
                    for i, line in enumerate (myfile):
                        if i>3:
                            data.append(line.rstrip("\r\n").split("\t"))
                    for i, (line ,line1,line2,line3) in \
                        enumerate(zip(data[0],data[1],data[2],data[3])):
                        total_read=line.split(" ")[2]
                        too_short=line1.split(" ")[4]
                        trim1=line2.split(" ")[1]
                        trim2=line3.split(" ")[1]
   
    #Sample page  
        page="\\begin{landscape} \n\
\subsection{All Samples} \n\
\\begin{minipage}[c][6in]{9.3in} \n\
\centering \n\
\\begin{minipage}[c][6in]{4in}\n\
\centering\n\
\\begin{tabular}{lrl}\n\
\multicolumn{3}{l}{Clipping/Trimming} \\\ \n\
\n\
\hline \n\
Total read              & "+ total_read +" \\\ \n\
Too short after clip    & "+ too_short +"  \\\ \n\
Trimmed R1              & "+ trim1 +" \\\ \n\
Trimmed R2              & "+ trim2 +"  \\\ \n\
\hline \n\
\\\\\\\\ \n\
\multicolumn{3}{l}{HS Metrics} \\\ \n\
\hline \n\
Total Reads             & "+Total_Reads+"  \\\ \n\
Total Duplicate         & "+Total_Duplicate+"\%  \\\ \n\
Duplicate Percentage    & "+Duplicate_Percentage+"\%   \\\ \n\
Optical Duplicate       & "+Optical_Duplicate+"\%   \\\ \n\
PCR Duplicate           & "+PCR_Duplicate+"\%   \\\ \n\
Off-Bait                & "+Off_Bait+"\%   \\\ \n\
Mean Bait Coverage      & "+Mean_Bait_Coverage+"   \\\ \n\
Mean Target Coverage    & "+Mean_Target_Coverage+"   \\\ \n\
Fold Enrichement        & "+Fold_Enrichement+"   \\\ \n\
Zero CVG Target         & "+Zero_CVG_Target+"\%  \\\ \n\
PCT-TARGET-BASES-2X     & "+PCT_TARGET_BASES_2X+"\%   \\\ \n\
PCT-TARGET-BASES-10X    & "+PCT_TARGET_BASES_10X+"\%   \\\ \n\
PCT-TARGET-BASES-20X    & "+PCT_TARGET_BASES_20X+"\%   \\\ \n\
PCT-TARGET-BASES-30X    & "+PCT_TARGET_BASES_30X+"\%   \\\ \n\
PCT-TARGET-BASES-40X    & "+PCT_TARGET_BASES_40X+"\%   \\\ \n\
PCT-TARGET-BASES-50X    & "+PCT_TARGET_BASES_50X+"\%   \\\ \n\
PCT-TARGET-BASES-100X   & "+PCT_TARGET_BASES_100X+"\%   \\\ \n\
\hline \n\
\\\\ \n\
\hline \n\
Mean Insert Size        &        & $\pm$ 10 \\\ \n\
Median Insert Size      &        & $\pm$ 9 \\\ \n\
Insert Size Std Dev     &         & $\pm$ 3 \\\ \n\
\hline \n\
\end{tabular}\n\
\end{minipage}%\n\
\\begin{minipage}[c][6in]{5in}\n\
\centering \n\
\includegraphics[width=4.5in]{"+coverage_graph_multi+"}\\\ \n\
\\vfill \n\
\includegraphics[width=4.5in]{"+coverage_graph_multi+"} \n\
\end{minipage} \n\
\end{minipage} \n\
\end{landscape} \n\
\n \
"
    #    print (page)       
        
        pages.append(page) 
    
    
    
    start="\documentclass[10pt,twoside,english]{scrartcl}\n\
\n\
\\usepackage[T1]{fontenc} \n\
\\usepackage[utf8]{inputenc} \n\
\n\
\\usepackage[letterpaper]{geometry} \n\
\geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm,footskip=1cm}\
\n\
\n\
\\usepackage{fancyhdr} \n\
\pagestyle{fancy} \n\
\n\
\setlength{\parindent}{0cm} \n\
\n\
\\usepackage{color} \n\
\n\
\\usepackage{babel} \n\
\n\
\\usepackage[unicode=true, \n \
bookmarks=true,bookmarksnumbered=false,bookmarksopen=false, \n \
breaklinks=false,pdfborder={0 0 0},backref=false,colorlinks=true] \n \
{hyperref} \n\
\hypersetup{pdftitle={Automatic Sequencing Report}, \n\
pdfauthor={Louis-Philippe Lemieux Perreault}, \n\
linkcolor={link_blue},citecolor={link_blue},urlcolor={link_blue}} \n\
\n\
\makeatletter \n\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User specified LaTeX commands. \n\
\\usepackage[activate={true,nocompatibility},final,tracking=true,kerning=true,spacing=true,factor=1100,stretch=10,shrink=10]\
{microtype} \n\
% activate={true,nocompatibility} - activate protrusion and expansion \n\
% final - enable microtype; use \"draft\" to disable \n\
% tracking=true, kerning=true, spacing=true - activate these techniques \n\
% factor=1100 - add 10% to the protrusion amount (default is 1000) \n\
% stretch=10, shrink=10 - reduce stretchability/shrinkability (default is\
20/20)\n\
\n\
% Some color definitions \n\
\\usepackage{xcolor} \n\
\definecolor{blue}{HTML}{0099CC} \n\
\definecolor{red}{HTML}{CC0000} \n\
\definecolor{green}{HTML}{669900} \n\
\definecolor{link_blue}{HTML}{21759B} \n\
\n\
% Caption and subcaptions \n\
\\usepackage[bf]{caption} \n\
\\usepackage{subcaption} \n\
\n\
% We want images \n\
\\usepackage{float} \n\
\\usepackage{graphicx} \n\
\n \
% Enables page in landscape mode \n\
\\usepackage{pdflscape} \n\
\n \
% Color in tables + dashed line \n\
\\usepackage{colortbl} \n\
\\usepackage{arydshln} \n\
\n \
% The Report number... \n\
\\newcommand{\RunName}{140430\_SN1064\_0137\_AC3TMHACXX} \n\
\n\
% Fancy header \n\
\setlength{\headheight}{18pt} \n\
\pagestyle{fancy} \n\
\\fancyhf{} \n\
\\fancyhead[RO,LE]{\includegraphics\
[height=0.78\headheight]{images/pgx_logo_small.png}}\n\
\\fancyfoot[RO,LE]{{\\footnotesize Automatic Sequencing Report}}\n\
\\fancyfoot[LO,RE]{{\\footnotesize\RunName}}\n\
\\fancyfoot[CO,CE]{{\\footnotesize\\thepage}}\n\
\\renewcommand{\headrulewidth}{1pt}\n\
\\renewcommand{\\footrulewidth}{1pt}\n\
\n\
% The plain fancy header \n\
\\fancypagestyle{plain}{\n\
    \\fancyhf{} % clear all header and footer fields\n\
        \\renewcommand{\headrulewidth}{0pt}\n\
            \\renewcommand{\\footrulewidth}{0pt}\n\
}\n\
\n\
% Tables numerotations are in Roman numbers \n\
\\renewcommand{\\thetable}{\Roman{table}}\n\
\n\
% Subfigure numbering\n\
\\renewcommand{\\thesubfigure}{\Alph{subfigure}}\n\
\n\
\let\\tempone\itemize\n\
\let\\temptwo\enditemize\n\
\\renewenvironment{itemize}{\\tempone\setlength{\itemsep}{0pt}}{\\temptwo}\n\
\n\
\let\\ttempone\enumerate\n\
\let\\ttemptwo\endenumerate\n\
\\renewenvironment{enumerate}{\\ttempone\setlength{\itemsep}{0pt}}{\\ttemptwo}\n\
\n\
\let\\tttempone\description\n\
\let\\tttemptwo\enddescription\n\
\\renewenvironment{description}{\\tttempone\setlength{\itemsep}{0pt}}{\\tttemptwo}\n\
\n\
\makeatother\n\
\n\
\\begin{document}\n\
\n\
% The title page\n\
\\title{Sequencing Report}\n\
\subtitle{\RunName}\n\
\\author{\includegraphics[width=0.6\\textwidth]{images/pgx_logo.pdf}}\n\
\date{\\today}\n\
\maketitle\n\
\n\
% The table of contents\n\
\microtypesetup{protrusion=false} % disables protrusion locally in the\
document\n\
\\tableofcontents{}\n\
\listoffigures{}\n\
\microtypesetup{protrusion=true} % enables protrusion\n\
\cleardoublepage\n\
\\newpage\n\
\n\
\section{Overview}\n\
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Mauris sed leo dolor.\n\
Curabitur sit amet tortor tempor orci sagittis interdum sit amet a odio.\n\
dapibus orci et convallis molestie. Praesent volutpat risus sed nibh elementum\n\
ullamcorper. Morbi vitae erat nibh. Integer malesuada urna non vestibulum\n\
sodales. Duis ornare dui urna, sed fermentum lectus faucibus elementum. Donec\n\
dui erat, dignissim sit amet neque nec, pretium tincidunt quam. Curabitur \n\
egestas tempus orci eu volutpat. Aliquam erat volutpat. Suspendisse ut ipsu \n\
accumsan, venenatis neque eget, porta ante. Figure~\\ref{pipeline_overview} \n\
the schematic of the pipeline used for run \RunName. \n\
\n\
\\begin{figure}[H]\n\
\centering \n\
\includegraphics[width=0.6\\textwidth]{images/flowchart.png} \n\
\caption[Pipeline overview]{Overview of the pipeline used for run \n\
\RunName.\label{pipeline_overview}} \n\
\end{figure} \n\
\n\
"

    end="% The summary for all samples \n\
\end{document} " 

    


#generate tex
    report = open("report.tex", "w")
    report.write(start)
    for  i,p in enumerate (pages):
        report.write (p)
    report.write (end )
# Close opend file
    report.close()
# generate pdf
    os.system("rubber -fds report.tex") 
# remove intermediate file
    for i in  (rmlist):
        os.remove(i)


# The parser object
desc = "Generate automatic report (version {}).".format(prog_version)
parser = argparse.ArgumentParser(description=desc)

# The input file
parser.add_argument("-v", "--version", action="version",
                    version="%(prog)s version {}".format(prog_version))

group = parser.add_argument_group("Pipeline Configuration")
group.add_argument("-i", "--input", type=str, metavar="FILE",
                   default="input_files.txt",
                   help=("A file containing the pipeline input files (one "
                         "sample per line, one or more file per sample "
                         "[%(default)s]"))
group.add_argument("-p", "--pipeline-config", type=str, metavar="FILE",
                   default="pipeline.conf",
                   help="The pipeline configuration file [%(default)s]")


try:
    # Getting and checking the options
    args = parse_args()
    check_args(args)

    # Summarize the options used
    print()
    print("{} \\".format(sys.argv[0]))
    for option, value in vars(args).items():
        option = option.replace("_", "-")
        print("    --{} {} \\".format(option, value))
    print()

    # Checking the input files
    input_files,sample_list = check_input_files(args.input)

    #Print report
    # Getting the pipeline steps
    what_to_run = get_pipeline_steps(args.pipeline_config)
    print_report (what_to_run,sample_list)

except KeyboardInterrupt:
    print >>sys.stderr, "Cancelled by user"
    sys.exit(0)
except ProgramError as e:
    parser.error(e.message)

