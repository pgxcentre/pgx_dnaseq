#!/usr/bin/env python3

import os
import re
import sys
import glob
import shutil
import logging
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pgx_dna_seq import ProgramError
from pgx_dna_seq.read_config import get_pipeline_steps


prog_version = "0.1"


def main():
    """The main function."""
    # Creating the option parser
    desc = "Generate automatic report (version {}).".format(prog_version)
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
        print_report(sample_list, steps)

    except KeyboardInterrupt:
        logging.info("Cancelled by user")
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    except ProgramError as e:
        logging.error(e)
        parser.error(e.message)

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
    # The required regular expression
    split_re = re.compile(r"\s+")
    sample_re = re.compile(r"^\w+")

    # Reading the input file
    all_samples = set()
    with open(filename, "r") as i_file:
        for line in i_file:
            samples = split_re.split(line.rstrip("\n"))
            samples = [os.path.basename(i) for i in samples]

            for sample in samples:
                if sample.endswith(".fastq") or sample.endswith(".fastq.gz"):
                    all_samples.add(re.search(r"(^\w+)_R[12]", sample).group(1))

                else:
                    all_samples.add(re.search(r"\w+", sample).group())

    logging.info("Found {:,d} samples".format(len(all_samples)))
    return sorted(all_samples)


def print_report(sample_list, pipeline_steps):
    """Creates the report."""
    pass


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
   
    return True


def parse_args(parser):
    """Parses the command line options and arguments."""
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s version {}".format(prog_version))
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

    return parser.parse_args()


if __name__ == "__main__":
    main()

## 
## 
## 
## 
## def std_dev(list):
##     [float(i) for i in list]
##     m=sum(list)/len(list)
##     z=sum(([(x-m)**2 for x in list]))/len(list)
##     return z**0.5
## 
## 
## def parse_args():
##     """Parses the command line options and arguments.
## 
##     :returns: A :py:class:`argparse.Namespace` object created by the
##               :py:mod:`argparse` module. It contains the values of the different
##               options.
## 
##     ===============   =======  ================================================
##         Options        Type                      Description
##     ===============   =======  ================================================
##     ===============   =======  ================================================
## 
##     .. note::
##         No option check is done here (except for the one automatically done by
##         :py:mod:`argparse`). Those need to be done elsewhere (see
##         :py:func:`checkArgs`).
## 
##     """
##     return parser.parse_args()
## 
## 
## def print_report(what_to_run,sample_list):
##     isize_path=""
##     isize_list =list ()
##     rmlist =['report.toc','report.lof','report.out','report.log','report.aux']
##     pages=list ()
##     val_total_read= list ()
##     val_too_short= list ()
##     val_trim1= list ()
##     val_trim2= list ()
##     val_Total_Reads= list ()
##     val_Total_Duplicate= list ()
##     val_Duplicate_Per= list ()
##     val_Optical_Duplicate= list ()
##     val_PCR_Duplicate= list ()
##     val_Off_Bait= list ()
##     val_Mean_Bait_Coverage= list ()
##     val_Mean_T_Coverage= list ()
##     val_Fold_Enrichement= list ()
##     val_Zero_CVG_Target= list ()
##     val_PCT_T_BASES_2X= list ()
##     val_PCT_T_BASES_10X= list ()
##     val_PCT_T_BASES_20X= list ()
##     val_PCT_T_BASES_30X= list ()
##     val_PCT_T_BASES_40X= list ()
##     val_PCT_T_BASES_50X= list ()
##     val_PCT_T_BASES_100X= list ()
##     val_Mean_Insert_Size= list ()
##     val_Median_Insert_Size= list ()
##     val_Standard_Deviation =list ()
##     for sample in  (sample_list):   
## 
##         for job_index, (job, job_options) in enumerate(what_to_run):
##             output_dir = os.path.join("output","{:02d}_{}".format(job_index +1,job.get_tool_name()))
## 
##         #Hsmetrics
##             if job.get_tool_name() ==  'HsMetrics' :
##                     hsmetrics_file=glob.glob(output_dir +'/'+sample+'*.hsmetrics')[0]
##                     myfile=open (hsmetrics_file)
##                     data= list ()
##                     for i, line in enumerate (myfile):
##                         if i==7:
##                             data.append(line.replace("_","-").rstrip("\r\n").split("\t"))
##                             Total_Reads=data[0][5]
##                             Off_Bait=data[0][15]                
##                             Mean_Bait_Coverage=data[0][20]      
##                             Fold_Enrichement=data[0][24]        
##                             Mean_T_Coverage=data[0][21]    
##                             Zero_CVG_Target=data[0][25]         
##                             PCT_T_BASES_2X=data[0][27]     
##                             PCT_T_BASES_10X=data[0][28]    
##                             PCT_T_BASES_20X=data[0][29]    
##                             PCT_T_BASES_30X=data[0][30]    
##                             PCT_T_BASES_40X=data[0][31]    
##                             PCT_T_BASES_50X=data[0][32]    
##                             PCT_T_BASES_100X=data[0][33]   
##                             
##                             val_Total_Reads.append(float(Total_Reads))
##                             val_Off_Bait.append(float(Off_Bait))
##                             val_Mean_Bait_Coverage.append(float(Mean_Bait_Coverage))
##                             val_Fold_Enrichement.append(float(Fold_Enrichement))
##                             val_Mean_T_Coverage.append(float(Mean_T_Coverage))
##                             val_Zero_CVG_Target.append(float(Zero_CVG_Target))
##                             val_PCT_T_BASES_2X.append(float(PCT_T_BASES_2X))
##                             val_PCT_T_BASES_10X.append(float(PCT_T_BASES_10X))
##                             val_PCT_T_BASES_20X.append(float(PCT_T_BASES_20X))
##                             val_PCT_T_BASES_30X.append(float(PCT_T_BASES_30X))
##                             val_PCT_T_BASES_40X.append(float(PCT_T_BASES_40X))
##                             val_PCT_T_BASES_50X.append(float(PCT_T_BASES_50X))
##                             val_PCT_T_BASES_100X.append(float(PCT_T_BASES_100X))
## 
##          #InserSize
##             if job.get_tool_name() ==  'InsertSize' :
##                 isize_file=glob.glob(output_dir +'/'+sample+'*.insertsize')[0]
##                 isize_graph=glob.glob(output_dir +'/'+sample+'*.png')[0]
##                 isize_path=output_dir +'/all_samples.png'
##                 myfile=open (isize_file)
##                 data= list ()
##                 for i, line in enumerate (myfile):
##                     if i==7:
##                         data.append(line.replace("_","-").rstrip("\r\n").split("\t"))
##                         Mean_Insert_Size=data[0][4]
##                         Median_Insert_Size=data[0][0]
##                         Standard_Deviation=data[0][5]
##                             
##                         val_Mean_Insert_Size.append(float(data[0][4]))
##                         val_Median_Insert_Size.append(float(data[0][0]))
##                         val_Standard_Deviation.append(float(data[0][5]))
##         
##         #Duplicate 
##             if job.get_tool_name() ==  'MarkDuplicates' :
##                 dedup_file=glob.glob(output_dir +'/'+sample+'*.dedup')[0]
##                 myfile=open (dedup_file)
##                 data= list ()
##                 for i, line in enumerate (myfile):
##                     if i==7:
##                         data.append(line.replace("_","-").rstrip("\r\n").split("\t"))
##                         Total_Duplicate=data[0][5]
##                         Optical_Duplicate=str (int(data[0][6]) * 100 / float(Total_Duplicate))
##                         PCR_Duplicate=str(100-float(Optical_Duplicate))
##                         Duplicate_Per=str(float(data[0][7]) * 100 )
## 
##                         val_Total_Duplicate.append(float(data[0][5]))
##                         val_Optical_Duplicate.append(int(data[0][6])*100 / int(Total_Duplicate))
##                         val_PCR_Duplicate.append(100-float(Optical_Duplicate))
##                         val_Duplicate_Per.append(float(data[0][7]) * 100)
##         
##        #Coverage_graph
##             if job.get_tool_name() ==  'CoverageGraph' :
##                     coverage_graph = glob.glob(output_dir +'/'+sample+'*.png')[0]
##                     coverage_graph=coverage_graph.replace(".png",".png").replace("/"+sample,"/"+sample)
##                     if not os.path.isfile(output_dir +'/'+sample+'.png'):
##                         shutil.copy(coverage_graph, output_dir +'/'+sample+'.png')
##                     coverage_graph=output_dir +'/'+sample+'.png'
##                     rmlist.append(coverage_graph)
##         #CoverageGraph_Multi
##             if (job.get_tool_name() ==  'CoverageGraph_Multi' and 'coverage_graph_multi'  not  in locals()):
##                     coverage_graph_multi = glob.glob(output_dir+'/all_samples*.png')[0]
##                     if not os.path.isfile(output_dir +'/all_samples.png'):
##                         shutil.copy(coverage_graph_multi, output_dir+'/all_samples.png')
##                     coverage_graph_multi=output_dir +'/all_samples.png'
##                     rmlist.append(coverage_graph_multi)
##          #ClipTrim 
##             if job.get_tool_name() ==  'ClipTrim' :
##                     clip_file=output_dir +'/'+sample+'.out'
##                     myfile=open (clip_file)
##                     data= list ()
##                     for i, line in enumerate (myfile):
##                         if i>3:
##                             data.append(line.rstrip("\r\n").split("\t"))
##                     for i, (line ,line1,line2,line3) in \
##                         enumerate(zip(data[0],data[1],data[2],data[3])):
##                         total_read=line.split(" ")[2]
##                         too_short=line1.split(" ")[4]
##                         trim1=line2.split(" ")[1]
##                         trim2=line3.split(" ")[1]
##                         
##                         val_total_read.append(int(line.split(" ")[2]))
##                         val_too_short.append(int(line1.split(" ")[4]))
##                         val_trim1.append(int(line2.split(" ")[1]))
##                         val_trim2.append(int(line3.split(" ")[1]))
##     #Sample page  
##         page="\\begin{landscape} \n\
## \subsection{"+sample+"} \n\
## \\begin{minipage}[c][6in]{9.3in} \n\
## \centering \n\
## \\begin{minipage}[c][6in]{4in}\n\
## \centering\n\
## \\begin{tabular}{lrl}\n\
## \multicolumn{3}{l}{Clipping/Trimming} \\\ \n\
## \n\
## \hline \n\
## Total read              & "+ total_read +" \\\ \n\
## Too short after clip    & "+ too_short +"  \\\ \n\
## Trimmed R1              & "+ trim1 +" \\\ \n\
## Trimmed R2              & "+ trim2 +"  \\\ \n\
## \hline \n\
## \\\\\\\\ \n\
## \multicolumn{3}{l}{HS Metrics} \\\ \n\
## \hline \n\
## Total Reads             & "+str("%.2f" % float(Total_Reads))+"  \\\ \n\
## Total Duplicate         & "+str("%.2f" % float( Total_Duplicate))+"  \\\ \n\
## Duplicate Percentage    & "+str("%.2f" % float(Duplicate_Per))+"\%   \\\ \n\
## Optical Duplicate       & "+str("%.2f" % float(Optical_Duplicate))+"\%   \\\ \n\
## PCR Duplicate           & "+str("%.2f" % float(PCR_Duplicate))+"\%   \\\ \n\
## Off-Bait                & "+str("%.2f" % float(Off_Bait))+"   \\\ \n\
## Mean Bait Coverage      & "+str("%.2f" % float(Mean_Bait_Coverage))+"   \\\ \n\
## Mean Target Coverage    & "+str("%.2f" % float(Mean_T_Coverage))+"   \\\ \n\
## Fold Enrichement        & "+str("%.2f" % float(Fold_Enrichement))+"   \\\ \n\
## Zero CVG Target         & "+str("%.2f" % float(Zero_CVG_Target))+"\%  \\\ \n\
## PCT-TARGET-BASES-2X     & "+str("%.2f" % float(PCT_T_BASES_2X))+"\%   \\\ \n\
## PCT-TARGET-BASES-10X    & "+str("%.2f" % float(PCT_T_BASES_10X))+"\%   \\\ \n\
## PCT-TARGET-BASES-20X    & "+str("%.2f" % float(PCT_T_BASES_20X))+"\%   \\\ \n\
## PCT-TARGET-BASES-30X    & "+str("%.2f" % float(PCT_T_BASES_30X))+"\%   \\\ \n\
## PCT-TARGET-BASES-40X    & "+str("%.2f" % float(PCT_T_BASES_40X))+"\%   \\\ \n\
## PCT-TARGET-BASES-50X    & "+str("%.2f" % float(PCT_T_BASES_50X))+"\%   \\\ \n\
## PCT-TARGET-BASES-100X   & "+str("%.2f" % float(PCT_T_BASES_100X))+"\%   \\\ \n\
## \hline \n\
## \\\\ \n\
## \hline \n\
## Mean Insert Size        & "+str("%.2f" % float(Mean_Insert_Size))+" \\\ \n\
## Median Insert Size      & "+str("%.2f" % float(Median_Insert_Size))+" \\\ \n\
## Insert Size Std Dev     & "+str("%.2f" % float(Standard_Deviation))+" \\\ \n\
## \hline \n\
## \end{tabular}\n\
## \end{minipage}%\n\
## \\begin{minipage}[c][6in]{5in}\n\
## \centering \n\
## \includegraphics[width=4.5in]{"+isize_graph+"}\\\ \n\
## \\vfill \n\
## \includegraphics[width=4.5in]{"+coverage_graph+"} \n\
## \end{minipage} \n\
## \end{minipage} \n\
## \end{landscape} \n\
## \n \
## "
##         
##         pages.append(page) 
##     
##     
##     
##     start="\documentclass[10pt,twoside,english]{scrartcl}\n\
## \n\
## \\usepackage[T1]{fontenc} \n\
## \\usepackage[utf8]{inputenc} \n\
## \n\
## \\usepackage[letterpaper]{geometry} \n\
## \geometry{verbose,tmargin=2cm,bmargin=2cm,lmargin=2cm,rmargin=2cm,footskip=1cm}\
## \n\
## \n\
## \\usepackage{fancyhdr} \n\
## \pagestyle{fancy} \n\
## \n\
## \setlength{\parindent}{0cm} \n\
## \n\
## \\usepackage{color} \n\
## \n\
## \\usepackage{babel} \n\
## \n\
## \\usepackage[unicode=true, \n \
## bookmarks=true,bookmarksnumbered=false,bookmarksopen=false, \n \
## breaklinks=false,pdfborder={0 0 0},backref=false,colorlinks=true] \n \
## {hyperref} \n\
## \hypersetup{pdftitle={Automatic Sequencing Report}, \n\
## pdfauthor={Automatic Sequencing Reporter}, \n\
## linkcolor={link_blue},citecolor={link_blue},urlcolor={link_blue}} \n\
## \n\
## \makeatletter \n\
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User specified LaTeX commands. \n\
## \\usepackage[activate={true,nocompatibility},final,tracking=true,kerning=true,spacing=true,factor=1100,stretch=10,shrink=10]\
## {microtype} \n\
## % activate={true,nocompatibility} - activate protrusion and expansion \n\
## % final - enable microtype; use \"draft\" to disable \n\
## % tracking=true, kerning=true, spacing=true - activate these techniques \n\
## % factor=1100 - add 10% to the protrusion amount (default is 1000) \n\
## % stretch=10, shrink=10 - reduce stretchability/shrinkability (default is\
## 20/20)\n\
## \n\
## % Some color definitions \n\
## \\usepackage{xcolor} \n\
## \definecolor{blue}{HTML}{0099CC} \n\
## \definecolor{red}{HTML}{CC0000} \n\
## \definecolor{green}{HTML}{669900} \n\
## \definecolor{link_blue}{HTML}{21759B} \n\
## \n\
## % Caption and subcaptions \n\
## \\usepackage[bf]{caption} \n\
## \\usepackage{subcaption} \n\
## \n\
## % We want images \n\
## \\usepackage{float} \n\
## \\usepackage{graphicx} \n\
## \n \
## % Enables page in landscape mode \n\
## \\usepackage{pdflscape} \n\
## \n \
## % Color in tables + dashed line \n\
## \\usepackage{colortbl} \n\
## \\usepackage{arydshln} \n\
## \n \
## % The Report number... \n\
## \\newcommand{\RunName}{"+run_name+"} \n\
## \n\
## % Fancy header \n\
## \setlength{\headheight}{18pt} \n\
## \pagestyle{fancy} \n\
## \\fancyhf{} \n\
## \\fancyhead[RO,LE]{\includegraphics\
## [height=0.78\headheight]{images/pgx_logo_small.png}}\n\
## \\fancyfoot[RO,LE]{{\\footnotesize Automatic Sequencing Report}}\n\
## \\fancyfoot[LO,RE]{{\\footnotesize\RunName}}\n\
## \\fancyfoot[CO,CE]{{\\footnotesize\\thepage}}\n\
## \\renewcommand{\headrulewidth}{1pt}\n\
## \\renewcommand{\\footrulewidth}{1pt}\n\
## \n\
## % The plain fancy header \n\
## \\fancypagestyle{plain}{\n\
##     \\fancyhf{} % clear all header and footer fields\n\
##         \\renewcommand{\headrulewidth}{0pt}\n\
##             \\renewcommand{\\footrulewidth}{0pt}\n\
## }\n\
## \n\
## % Tables numerotations are in Roman numbers \n\
## \\renewcommand{\\thetable}{\Roman{table}}\n\
## \n\
## % Subfigure numbering\n\
## \\renewcommand{\\thesubfigure}{\Alph{subfigure}}\n\
## \n\
## \let\\tempone\itemize\n\
## \let\\temptwo\enditemize\n\
## \\renewenvironment{itemize}{\\tempone\setlength{\itemsep}{0pt}}{\\temptwo}\n\
## \n\
## \let\\ttempone\enumerate\n\
## \let\\ttemptwo\endenumerate\n\
## \\renewenvironment{enumerate}{\\ttempone\setlength{\itemsep}{0pt}}{\\ttemptwo}\n\
## \n\
## \let\\tttempone\description\n\
## \let\\tttemptwo\enddescription\n\
## \\renewenvironment{description}{\\tttempone\setlength{\itemsep}{0pt}}{\\tttemptwo}\n\
## \n\
## \makeatother\n\
## \n\
## \\begin{document}\n\
## \n\
## % The title page\n\
## \\title{Sequencing Report}\n\
## \subtitle{\RunName}\n\
## \\author{\includegraphics[width=0.6\\textwidth]{images/pgx_logo.pdf}}\n\
## \date{\\today}\n\
## \maketitle\n\
## \n\
## % The table of contents\n\
## \microtypesetup{protrusion=false} % disables protrusion locally in the\
## document\n\
## \\tableofcontents{}\n\
## \listoffigures{}\n\
## \microtypesetup{protrusion=true} % enables protrusion\n\
## \cleardoublepage\n\
## \section{Overview}\n\
## This report includes metrics collected using picard Hsmetrics , picard\n\
## MarkDuplicate , picard insert size metrics . \n\
##  Figure~\\ref{pipeline_overview} \n\
## the schematic of the pipeline used for run \RunName. \n\
## \n\
## \\begin{figure}[H]\n\
## \centering \n\
## \includegraphics[width=0.6\\textwidth]{images/flowchart.pdf} \n\
## \caption[Pipeline overview]{Overview of the pipeline used for run \n\
## \RunName.\label{pipeline_overview}} \n\
## \end{figure} \n\
## \n\
## "
## 
##     end="% The summary for all samples \n\
## \end{document} " 
## 
##     sample_num=len(sample_list)
## 
## #Generate Insert Size plot
##     plt.xlabel('Sample')
##     plt.ylabel('Mean Insert Size')
##     plt.title('Histogram Mean Insert size values')
##     x=list (range(1,sample_num+1))
##     plt.axis([0,sample_num+1,0,500])
##     plt.plot(x,val_Mean_Insert_Size,color='red',marker='o', markersize=10)
##     plt.savefig(isize_path)
## 
## #all samples page
##     sample_num=len(sample_list)
##     all_samples="\\begin{landscape} \n\
## \subsection{all\_sample} \n\
## \\begin{minipage}[c][6in]{9.3in} \n\
## \centering \n\
## \\begin{minipage}[c][6in]{4in}\n\
## \centering\n\
## \\begin{tabular}{lrl}\n\
## \multicolumn{3}{l}{Clipping/Trimming} \\\ \n\
## \n\
## \hline \n\
## Total read              & "+ str("%.2f"  %(sum(val_total_read)/sample_num)) +"& $\pm$ "+str("%.2f"  % std_dev(val_total_read))+"  \\\ \n\
## Too short after clip    & "+ str("%.2f" % (sum(val_too_short)/sample_num)) +" & $\pm$ "+str("%.2f"  % std_dev(val_too_short))+"\\\ \n\
## Trimmed R1              & "+ str("%.2f" % (sum(val_trim1)/sample_num)) +" & $\pm$ "+str("%.2f"  % std_dev(val_trim1))+"\\\ \n\
## Trimmed R2              & "+ str("%.2f" % (sum(val_trim2)/sample_num)) +" & $\pm$ "+str("%.2f"  % std_dev(val_trim2))+" \\\ \n\
## \hline \n\
## \\\\\\\\ \n\
## \multicolumn{3}{l}{HS Metrics} \\\ \n\
## \hline \n\
## Total Reads             & "+str("%.2f" % (sum(val_Total_Reads)/sample_num))+" & $\pm$ "+str("%.2f"  % std_dev(val_Total_Reads))+"\\\ \n\
## Total Duplicate         & "+str("%.2f" %(sum(val_Total_Duplicate)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_Total_Duplicate))+"\\\ \n\
## Duplicate Percentage    & "+str("%.2f" %(sum(val_Duplicate_Per)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_Duplicate_Per))+" \\\ \n\
## Optical Duplicate       & "+str("%.2f" % (sum(val_Optical_Duplicate)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_Optical_Duplicate))+" \\\ \n\
## PCR Duplicate           & "+str("%.2f" % (sum(val_PCR_Duplicate)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_PCR_Duplicate))+" \\\ \n\
## Off-Bait                & "+str("%.2f" % (sum(val_Off_Bait)/sample_num))+"  & $\pm$ "+str("%.2f"  % std_dev(val_Off_Bait))+" \\\ \n\
## Mean Bait Coverage      & "+str("%.2f" % (sum(val_Mean_Bait_Coverage)/sample_num))+"  & $\pm$ "+str("%.2f"  % std_dev(val_Mean_Bait_Coverage))+" \\\ \n\
## Mean Target Coverage    & "+str("%.2f" % (sum(val_Mean_T_Coverage)/sample_num))+"  & $\pm$ "+str("%.2f"  % std_dev(val_Mean_T_Coverage))+" \\\ \n\
## Fold Enrichement        & "+str("%.2f" % (sum(val_Fold_Enrichement)/sample_num))+"  & $\pm$ "+str("%.2f"  % std_dev(val_Fold_Enrichement))+" \\\ \n\
## Zero CVG Target         & "+str("%.2f" % (sum(val_Zero_CVG_Target)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_Zero_CVG_Target))+" \\\ \n\
## PCT-TARGET-BASES-2X     & "+str("%.2f" % (sum(val_PCT_T_BASES_2X)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_PCT_T_BASES_2X))+" \\\ \n\
## PCT-TARGET-BASES-10X    & "+str("%.2f" % (sum(val_PCT_T_BASES_10X)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_PCT_T_BASES_10X))+" \\\ \n\
## PCT-TARGET-BASES-20X    & "+str("%.2f" % (sum(val_PCT_T_BASES_20X)/sample_num))+"\%  & $\pm$  "+str("%.2f"  % std_dev(val_PCT_T_BASES_20X))+"\\\ \n\
## PCT-TARGET-BASES-30X    & "+str("%.2f" %(sum(val_PCT_T_BASES_30X)/sample_num))+"\%  & $\pm$ "+str("%.2f"  %std_dev(val_PCT_T_BASES_30X))+" \\\ \n\
## PCT-TARGET-BASES-40X    & "+str("%.2f" %(sum(val_PCT_T_BASES_40X)/sample_num))+"\%  & $\pm$ "+str("%.2f"  % std_dev(val_PCT_T_BASES_40X))+" \\\ \n\
## PCT-TARGET-BASES-50X    & "+str("%.2f"%(sum(val_PCT_T_BASES_50X)/sample_num))+"\%  & $\pm$ "+str("%.2f"  %std_dev(val_PCT_T_BASES_50X))+" \\\ \n\
## PCT-TARGET-BASES-100X   & "+str("%.2f"%(sum(val_PCT_T_BASES_100X)/sample_num))+"\%  & $\pm$ "+str("%.2f" %std_dev(val_PCT_T_BASES_100X))+"\\\ \n\
## \hline \n\
## \\\\ \n\
## \hline \n\
## Mean Insert Size        & "+str("%.2f"%(sum(val_Mean_Insert_Size)/sample_num))+"  & $\pm$ "+str("%.2f" %std_dev(val_Mean_Insert_Size))+" \\\ \n\
## Median Insert Size      & "+str("%.2f"%(sum(val_Median_Insert_Size)/sample_num))+"  & $\pm$ "+str("%.2f" %std_dev(val_Median_Insert_Size))+" \\\ \n\
## Insert Size Std Dev     & "+str("%.2f"%(sum(val_Standard_Deviation)/sample_num))+"  & $\pm$ "+str("%.2f" % std_dev(val_Standard_Deviation))+" \\\ \n\
## \hline \n\
## \end{tabular}\n\
## \end{minipage}%\n\
## \\begin{minipage}[c][6in]{5in}\n\
## \centering \n\
## \includegraphics[width=4.5in]{"+isize_path+"}\\\ \n\
## \\vfill \n\
## \includegraphics[width=4.5in]{"+coverage_graph_multi+"} \n\
## \end{minipage} \n\
## \end{minipage} \n\
## \end{landscape} \n\
## \n \
## "
## #generate tex
##     report = open("report.tex", "w")
##     report.write(start)
##     report.write(all_samples)
##     for  i,p in enumerate (pages):
##         report.write (p)
##     report.write (end )
## # Close opend file
##     report.close()
## # generate pdf
##     os.system("rubber -fds report.tex") 
## # remove intermediate file
##     for i in  (rmlist):
##         os.remove(i)
## 
## 
## 
## try:
## 
## except KeyboardInterrupt:
##     print >>sys.stderr, "Cancelled by user"
##     sys.exit(0)
## except ProgramError as e:
##     parser.error(e.message)
