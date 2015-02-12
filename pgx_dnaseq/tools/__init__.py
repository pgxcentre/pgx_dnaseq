
# This file is part of pgx_dnaseq
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import os
import shlex
from glob import glob
from math import ceil
from tempfile import NamedTemporaryFile
from subprocess import check_call, SubprocessError

from .. import ProgramError


__all__ = ["bwa", "fastq_mcf", "fastqc", "gatk", "picard_tools", "samtools",
           "bowtie2", "bcftools", "pgx_coverage_graph",
           "pgx_read_quality_graph"]


class GenericTool(object):

    # The type that can be use for options
    INPUT = 1
    INPUTS = 2
    OUTPUT = 3
    OPTIONAL = 4
    REQUIREMENT = 5
    INPUT_TO_SPLIT = 6

    # By default, a tool produces usable data
    _produce_data = True

    # By default, no multiple inputs
    _merge_all_inputs = False

    # The local tool configuration
    _tool_configuration = {}

    # By default, we run locally
    _locally = True

    # The script preamble
    _script_preamble = ""

    def __init__(self):
        """Initialize an new GeneticTool object."""
        # The generic command for the generic tool
        pass

    def produce_usable_data(self):
        """Returns True if the tool produces data. False otherwise."""
        return self._produce_data

    def need_to_merge_all_inputs(self):
        """Returns True if the tool has multiple inputs. False otherwise."""
        return self._merge_all_inputs

    @staticmethod
    def set_tool_configuration(drmaa_options):
        """Sets the configuration for all the tools."""
        GenericTool._tool_configuration = drmaa_options

    @staticmethod
    def get_tool_configuration():
        """Get the configuration for all the tools."""
        return GenericTool._tool_configuration

    @staticmethod
    def do_not_run_locally():
        """Do not run the tools locally (sets _locally to False)."""
        GenericTool._locally = False

    @staticmethod
    def set_script_preamble(preamble):
        """Set the script preamble when using DRMAA."""
        GenericTool._script_preamble = preamble

    @staticmethod
    def get_script_preamble():
        """Returns the script preamble when using DRMAA."""
        return GenericTool._script_preamble

    @staticmethod
    def run_locally():
        """Do the tools need to be run locally or not."""
        return GenericTool._locally

    @staticmethod
    def get_tool_bin_dir(tool_name):
        """Returns the binary directory (empty string if none specified)."""
        # Getting all the tool configuration
        tool_conf = GenericTool.get_tool_configuration()

        # By default, the binary directory is an empty string (meaning that the
        # tool's binary directory is in the PATH variable)
        bin_dir = ""
        if (tool_name in tool_conf) and ("bin_dir" in tool_conf[tool_name]):
            bin_dir = tool_conf[tool_name]["bin_dir"]

        # Returning the binary directory
        return bin_dir

    def get_tool_name(self):
        """Returns the tool name."""
        try:
            return self._tool_name
        except AttributeError:
            return self.__class__.__name__

    def get_version(self):
        """Returns the version of the tool."""
        try:
            return self._version
        except AttributeError:
            m = "{}: no version".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_executable(self):
        """Returns the executable."""
        try:
            return self._exec
        except AttributeError:
            m = "{}: no executable".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_input_type(self):
        """Returns the input type of the tool."""
        try:
            return self._input_type
        except AttributeError:
            m = "{}: no input type".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_output_type(self):
        """Returns the output type of the tool."""
        try:
            return self._output_type
        except AttributeError:
            m = "{}: no output type".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_command(self):
        """Returns the command options."""
        try:
            return self._command
        except AttributeError:
            m = "{}: no command".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_stdout(self):
        """Returns the tool's STDOUT."""
        try:
            return self._stdout
        except AttributeError:
            m = "{}: no STDOUT".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_stderr(self):
        """Returns the tool's STDERR."""
        try:
            return self._stderr
        except AttributeError:
            m = "{}: no STDERR".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_required_options(self):
        """Returns the required options."""
        try:
            return self._required_options
        except AttributeError:
            m = "{}: required_options is None".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_suffix(self):
        """Returns the suffix."""
        try:
            return self._suffix
        except AttributeError:
            m = "{}: suffix is None".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def execute(self, tool_options, out_dir=None):
        """Executes the tool."""
        # The name of the job
        tool_name = self.get_tool_name()

        # Are we in need of a bulk submission?
        bulk, nb_chunks, file_to_split = GenericTool._is_bulk_job(
            GenericTool.get_tool_configuration(),
            tool_name,
        )

        # Do we need to split a file for bulk jobs?
        nb_split=None
        original_output_name = None
        if bulk:
            # Is there an output file?
            if "output" not in tool_options:
                m = "{}: cannot run in bulk job".format(tool_name)
                raise ProgramError(m)

            split_name, nb_split = GenericTool._split_file(
                file_to_split=tool_options[file_to_split],
                nb_chunks=nb_chunks,
                out_dir=out_dir,
                out_dir_suffix=tool_options["sample_id"],
            )

            # Changing the name of the split file
            tool_options[file_to_split] = split_name

            # There is now one output per spit job
            original_output_name = tool_options["output"]
            tool_options["output"] += "_$PBS_ARRAYID"

        # Checks the options
        checked_options = self.check_options(tool_options)

        # Create the command
        bin_dir = GenericTool.get_tool_bin_dir(tool_name)
        job_command = [os.path.join(bin_dir, self.get_executable())]
        job_command += self.get_command().format(**checked_options).split()

        # The STDOUT and STDERR files
        job_stdout = self.get_stdout().format(**checked_options)
        job_stderr = self.get_stderr().format(**checked_options)

        # Execute it
        if GenericTool.run_locally():
            GenericTool._execute_command_locally(
                command=job_command,
                stdout=job_stdout,
                stderr=job_stderr,
            )
        else:
            # Getting the tool walltime and nodes variable (for DRMAA)
            walltime, nodes = GenericTool._create_drmaa_var(
                GenericTool.get_tool_configuration(),
                tool_name,
            )

            if bulk:
                GenericTool._execute_bulk_command_drmaa(
                    preamble=GenericTool.get_script_preamble(),
                    command=job_command,
                    stdout=job_stdout,
                    stderr=job_stderr,
                    out_dir=out_dir,
                    job_name=tool_name,
                    walltime=walltime,
                    nodes=nodes,
                    nb_chunks=nb_split,
                )

                # Merging the bulk jobs
                try:
                    self.merge_bulk_results(original_output_name,
                                            tool_options["output"],
                                            nb_split)

                except AttributeError:
                    m = "{}: unable to join bulk results".format(tool_name)
                    raise ProgramError(m)

            else:
                GenericTool._execute_command_drmaa(
                    command=job_command,
                    stdout=job_stdout,
                    stderr=job_stderr,
                    out_dir=out_dir,
                    job_name=tool_name,
                    walltime=walltime,
                    nodes=nodes,
                    preamble=GenericTool.get_script_preamble(),
                )

    @staticmethod
    def _execute_command_locally(command, stdout=None, stderr=None):
        """Executes a command using the subprocess module."""
        # The stdout and stderr files
        if stdout is not None:
            stdout = open(stdout, "wb")
        if stderr is not None:
            stderr = open(stderr, "wb")

        # The process
        try:
            check_call(command, stdout=stdout, stderr=stderr)
        except SubprocessError:
            # Constructing the error message
            m = "The following command failed:\n\n"
            m += "    {}\n\n".format(" ".join(command))

            # The name of the log file
            log_filename = "log file"
            if stderr is not None:
                log_filename = stderr.name
            m += "Check {} for more detail".format(log_filename)

            # Raising the exception
            raise ProgramError(m)

        except FileNotFoundError:
            m = "{}: no such executable".format(command[0])
            raise ProgramError(m)

        finally:
            # Closing the output files
            if stdout is not None:
                stdout.close()
            if stderr is not None:
                stderr.close()

    @staticmethod
    def _execute_command_drmaa(preamble, command, stdout, stderr, out_dir,
                               job_name, walltime, nodes):
        """Executes a command using DRMAA."""
        # Creating the script in a temporary file
        tmp_file = NamedTemporaryFile(mode="w", suffix="_execute.sh",
                                      delete=False, dir=out_dir)

        # Writing the shebang
        print("#!/usr/bin/env bash", file=tmp_file)

        # Writing the preamble
        print(preamble, file=tmp_file)

        # Writing the command
        print(command[0], end=" ", file=tmp_file)
        for chunck in command[1:]:
            print(shlex.quote(chunck), end=" ", file=tmp_file)
        print("> {}".format(shlex.quote(stdout)), end=" ", file=tmp_file)
        print("2> {}".format(shlex.quote(stderr)), file=tmp_file, end="\n\n")

        # Closing the temporary file
        tmp_file.close()

        # Making the script executable
        os.chmod(tmp_file.name, 0o755)

        # Try executing the script using DRMAA
        try:
            import drmaa
        except ImportError:
            # Executing it locally
            GenericTool._execute_command_locally([tmp_file.name])
        else:
            # Initializing a new DRMAA session
            s = drmaa.Session()
            s.initialize()

            # Creating the job template
            job = s.createJobTemplate()
            job.remoteCommand = tmp_file.name
            job.jobName = "_{}".format(job_name)
            job.workingDirectory = os.getcwd()
            if walltime is not None:
                job.hardWallclockTimeLimit = walltime
            if nodes is not None:
                job.nativeSpecification = nodes

            # Running the job
            job_id = s.runJob(job)

            # Waiting for the job
            ret_val = s.wait(job_id, drmaa.Session.TIMEOUT_WAIT_FOREVER)

            # Deleting the job
            s.deleteJobTemplate(job)

            # Closing the connection
            s.exit()

            # Checking if there were problem
            if not GenericTool._is_job_completed(ret_val):
                m = "Could not run {}".format(tmp_file.name)
                raise ProgramError(m)

        # Removing the file
        os.remove(tmp_file.name)

    @staticmethod
    def _execute_bulk_command_drmaa(preamble, command, stdout, stderr, out_dir,
                                    job_name, walltime, nodes, nb_chunks):
        """Executes a bulk command using DRMAA."""
        # Creating the script in a temporary file
        tmp_file = NamedTemporaryFile(mode="w", suffix="_execute.sh",
                                      delete=False, dir=out_dir)

        # Writing the shebang
        print("#!/usr/bin/env bash", file=tmp_file)

        # Writing the preamble
        print(preamble, file=tmp_file)

        # Writing the command
        print(command[0], end=" ", file=tmp_file)
        for chunck in command[1:]:
            safe_chunk = shlex.quote(chunck)
            safe_chunk = safe_chunk.replace("$PBS_ARRAYID", "'${PBS_ARRAYID}'")
            print(safe_chunk, end=" ", file=tmp_file)
        print("> {}".format(shlex.quote(stdout)), end=" ", file=tmp_file)
        print("2> {}".format(shlex.quote(stderr)), file=tmp_file, end="\n\n")

        # Closing the temporary file
        tmp_file.close()

        # Making the script executable
        os.chmod(tmp_file.name, 0o755)

        # Try executing the script using DRMAA
        try:
            import drmaa
        except ImportError:
            # Executing it locally
            m = ("{} only work with DRMAA when using bulk "
                 "submission".format(job_name))
            raise ProgramError(m)

        else:
            # Initializing a new DRMAA session
            s = drmaa.Session()
            s.initialize()

            # Creating the job template
            job = s.createJobTemplate()
            job.remoteCommand = tmp_file.name
            job.jobName = "_{}".format(job_name)
            job.workingDirectory = os.getcwd()
            if walltime is not None:
                job.hardWallclockTimeLimit = walltime
            if nodes is not None:
                job.nativeSpecification = nodes

            # Running the job in array
            joblist = s.runBulkJobs(job, 1, nb_chunks, 1)

            # Waiting for the jobs
            s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, False)

            # Getting the return values
            ret_vals = [
                s.wait(j, drmaa.Session.TIMEOUT_WAIT_FOREVER) for j in joblist
            ]

            # Deleting the job
            s.deleteJobTemplate(job)

            # Closing the connection
            s.exit()

            # Checking if there were problem
            for ret_val in ret_vals:
                if not GenericTool._is_job_completed(ret_val):
                    m = "Could not run {}".format(tmp_file.name)
                    raise ProgramError(m)

        # Removing the file
        os.remove(tmp_file.name)

    @staticmethod
    def _split_file(file_to_split, nb_chunks, out_dir, out_dir_suffix):
        """Split a file to launch a bulk job."""
        # Now, we need to split the target file...
        to_split = None
        with open(file_to_split, "r") as i_file:
            to_split = i_file.read().splitlines()

        # Getting the number of lines per chunk
        nb_lines = ceil(len(to_split) / nb_chunks)

        # Just to be sure, we keep the number of files
        nb_files = 0

        # Getting the name and extension of the files
        name, ext = os.path.splitext(os.path.basename(file_to_split))

        # Splitting
        dirname = os.path.join(out_dir, "{}_chunks".format(out_dir_suffix))
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        filename = os.path.join(dirname, name + "_{i}" + ext)

        split_generator = (
            to_split[i:i+nb_lines] for i in range(0, len(to_split), nb_lines)
        )
        for i, lines in enumerate(split_generator):
            o_filename = filename.format(i=i+1)
            with open(o_filename, "w") as o_file:
                print(*lines, sep="\n", file=o_file)
            nb_files += 1

        # Returning the name of the split files
        return filename.format(i="$PBS_ARRAYID"), nb_files

    @staticmethod
    def _is_job_completed(job):
        """Checks the job status and return False if not completed."""
        if job.hasCoreDump or job.wasAborted or job.hasSignal:
            return False

        if job.exitStatus != 0:
            return False

        return True

    @staticmethod
    def _create_drmaa_var(options, job_name):
        """Creates "walltime" and "nodes" variables for the job."""
        # Creating the walltime variable
        walltime = None
        nodes = None
        if job_name in options:
            job_options = options[job_name]
            if "walltime" in job_options:
                walltime = bytes(job_options["walltime"], encoding="ascii")
            if ("nb_node" in job_options) and ("nb_proc" in job_options):
                nodes = "-l nodes={}:ppn={}".format(job_options["nb_node"],
                                                    job_options["nb_proc"])
                nodes = bytes(nodes, encoding="ascii")

        # Returns the walltime and the nodes information
        return walltime, nodes

    @staticmethod
    def _is_bulk_job(options, job_name):
        """Checks if the tool needs splitting (for array submission)."""
        bulk_submission = False
        nb_chunks = 1
        split_file = None
        if job_name in options:
            job_options = options[job_name]
            if ("nb_chunks" in job_options) and ("split_file" in job_options):
                nb_chunks = int(job_options["nb_chunks"])
                split_file = job_options["split_file"]
                bulk_submission = True

        return bulk_submission, nb_chunks, split_file

    def check_options(self, options):
        """Checks the tool options."""
        # Options that are automatically generated
        #   - input
        #   - output
        #   - prefix
        #   - sample_id

        # The "safe" options
        safe_options = {}

        # The required options
        required_options = self.get_required_options()

        # We check every required options
        for option_name, option_type in required_options.items():
            if (option_name not in options) and option_type != self.OPTIONAL:
                m = "{}: missing required option".format(option_name)
                raise ProgramError(m)

            # Checking if the option is an input file
            if option_type == self.INPUT:
                # The file should exists
                if not os.path.isfile(options[option_name]):
                    m = "{}: no such file".format(options[option_name])
                    raise ProgramError(m)

                # Option is now safe
                safe_options[option_name] = options[option_name]

            # Checking if the option is an input split file
            elif option_type == self.INPUT_TO_SPLIT:
                # Just checking that there are files
                globname = options[option_name]
                if len(glob(globname.replace("$PBS_ARRAYID", "*"))) < 1:
                    m = "{}: no such files".format(options[option_name])
                    raise ProgramError(m)

                # Option is now safe
                safe_options[option_name] = options[option_name]

            # Checking if the option is multiple input files
            elif option_type == self.INPUTS:
                # The files should all exists
                for filename in options[option_name]:
                    if not os.path.isfile(filename):
                        m = "{}: no such file".format(filename)
                        raise ProgramError(m)

                # They all exists, so we merge the input files using spaces
                safe_options[option_name] = " ".join(options[option_name])

            # Checking if the option is an output file
            elif option_type == self.OUTPUT:
                # We should be able to write the file
                test_filename = "{}.test".format(options[option_name])
                try:
                    with open(test_filename, 'w') as o_file:
                        pass
                except IOError:
                    m = "{}: cannot write file".format(options[option_name])
                    raise ProgramError(m)
                else:
                    if os.path.isfile(test_filename):
                        os.remove(test_filename)

                # Option is now safe
                safe_options[option_name] = options[option_name]

            # Checking the requirements
            elif option_type == self.REQUIREMENT:
                # There is nothing much to do here, since we don't know what
                # kind of option this is... It will be check by the tool on
                # runtime...
                safe_options[option_name] = options[option_name]

            # Checking if the option is optional
            elif option_type == self.OPTIONAL:
                # We only check if the option is available
                if option_name in options:
                    # Yes, so it *should* be safe
                    safe_options[option_name] = options[option_name]
                else:
                    # No, so we add an empty string to the safe options
                    safe_options[option_name] = ""

            # This option type is not valid...
            else:
                m = ("{}: {}: {}: not a valid option "
                     "type".format(self.__class__.__name__, option_name,
                                   option_type))
                raise ProgramError(m)

        # Returning the "safe" options
        return safe_options

    def read_report(self, *kargs, **kwargs):
        """Reads the report from one prefix."""
        raise NotImplementedError("No reporting options for "
                                  "'{}'".format(self.get_tool_name()))
