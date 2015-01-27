
# This file is part of pgx_dnaseq
#
# This work is licensed under The MIT License (MIT). To view a copy of this
# license, visit http://opensource.org/licenses/MIT


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = ("Copyright 2015 Beaulieu-Saucier Universite de Montreal "
                 "Pharmacogenomics Centre. All rights reserved.")
__license__ = "MIT"


import os
import shlex
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

    # By default, a tool produces usable data
    _produce_data = True

    # By default, no multiple inputs
    _merge_all_inputs = False

    # The local tool configuration
    __tool_configuration = {}

    # By default, we run locally
    __locally = True

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
        GenericTool.__tool_configuration = drmaa_options

    @staticmethod
    def get_tool_configuration():
        """Get the configuration for all the tools."""
        return GenericTool.__tool_configuration

    @staticmethod
    def do_not_run_locally():
        """Do not run the tools locally (sets __locally to False)."""
        GenericTool.__locally = False

    @staticmethod
    def run_locally():
        """Do the tools need to be run locally or not."""
        return GenericTool.__locally

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
            GenericTool.__execute_command_locally(job_command, job_stdout,
                                                  job_stderr)
        else:
            # Getting the tool walltime and nodes variable (for DRMAA)
            walltime, nodes = GenericTool.__create_drmaa_var(
                GenericTool.get_tool_configuration(),
                tool_name,
            )
            GenericTool.__execute_command_drmaa(job_command, job_stdout,
                                                job_stderr, out_dir, tool_name,
                                                walltime, nodes)

    @staticmethod
    def __execute_command_locally(command, stdout=None, stderr=None):
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
    def __execute_command_drmaa(command, stdout, stderr, out_dir, job_name,
                                walltime, nodes):
        """Executes a command using DRMAA."""
        # Creating the script in a temporary file
        tmp_file = NamedTemporaryFile(mode="w", suffix="_execute.sh",
                                      delete=False, dir=out_dir)
        print("#!/usr/bin/env bash", file=tmp_file)
        print(command[0], end=" ", file=tmp_file)
        for chunck in command[1:]:
            print(shlex.quote(chunck), end=" ", file=tmp_file)
        print("> {}".format(shlex.quote(stdout)), end=" ", file=tmp_file)
        print("2> {}".format(shlex.quote(stderr)), file=tmp_file)
        tmp_file.close()

        # Making the script executable
        os.chmod(tmp_file.name, 0o755)

        # Try executing the script using DRMAA
        try:
            import drmaa
        except ImportError:
            # Executing it locally
            GenericTool.__execute_command_locally([tmp_file.name])
        else:
            # Initializing a new DRMAA session
            s = drmaa.Session()
            s.initialize()

            # Creating the job template
            job = s.createJobTemplate()
            job.remoteCommand = tmp_file.name
            job.jobName = "_{}".format(job_name)
            job.workingDirectory = os.getcwd()
            job.jobEnvironment = os.environ
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
            if ret_val.exitStatus != 0:
                m = "Could not run {}".format(tmp_file.name)
                raise ProgramError(m)

        # Removing the file
        os.remove(tmp_file.name)

    @staticmethod
    def __create_drmaa_var(options, job_name):
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
