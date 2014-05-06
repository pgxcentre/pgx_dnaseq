__all__ = ["Java", "JAR"]

import os

from pgx_dna_seq import ProgramError
from pgx_dna_seq.tool import GenericTool


class Java(GenericTool):

    # The version of the tool
    _version = "1.7.0_51"

    # The executable
    _exec = "java"

    def __init__(self):
        """Initialize a Samtools instance."""
        pass


class JAR(Java):

    # The command specific to a jar
    _jar_command = "-Xmx{java_memory} {java_other_opt} -jar {jar_file}"

    # The jar required options
    # The description of the required options
    _jar_required_options = {"java_memory":    GenericTool.REQUIREMENT,
                             "jar_file":       GenericTool.REQUIREMENT,
                             "java_other_opt": GenericTool.OPTIONAL}

    def __init__(self):
        """Initialize a _SAM2BAM instance."""
        pass

    @staticmethod
    def get_tool_jar_dir(tool_name):
        """Returns the JAR directory (empty string if none specified)."""
        # Getting all the tool configuration
        tool_conf = GenericTool.get_tool_configuration()

        # By default, the JAR directory is None (meaning that the default
        # tool's JAR directory should be used)
        jar_dir = None
        if (tool_name in tool_conf) and ("jar_dir" in tool_conf[tool_name]):
            jar_dir = tool_conf[tool_name]["jar_dir"]

        # Returning the JAR directory
        return jar_dir

    def get_jar_file(self):
        """Returns the tool's jar file."""
        # Getting the optional jar directory
        jar_location = JAR.get_tool_jar_dir(self.get_tool_name())
        if jar_location is None:
            jar_location = self.get_jar_location()

        try:
            return os.path.join(jar_location, self._jar)
        except AttributeError:
            m = "{}: no jar".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_jar_required_options(self):
        """Returns the jar required options."""
        return self._jar_required_options

    def get_jar_location(self):
        """Gets the directory where the JAR is located."""
        try:
            return self._jar_location
        except AttributeError:
            m = ("{}: default jar location was not "
                 "set".format(self.__class__.__name__))
            raise ProgramError(m)

    def get_required_options(self):
        """Returns the required options."""
        try:
            return self.required_options
        except AttributeError:
            m = ("{}: required options were no "
                 "set".format(self.__class__.__name__))
            raise ProgramError(m)

    def get_child_required_options(self):
        """Returns the required options of the child."""
        try:
            return self._required_options
        except AttributeError:
            m = "{}: required_options is None".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def get_jar_command(self):
        """Returns the JAR command."""
        return self._jar_command

    def get_command(self):
        """Returns the command options."""
        try:
            return self.command
        except AttributeError:
            m = "{}: command was not set".format(self.__class__.__name__)
            raise ProgramError(m)

    def get_child_command(self):
        """Returns the command options of the child."""
        try:
            return self._command
        except AttributeError:
            m = "{}: no command".format(self.__class__.__name__)
            raise NotImplementedError(m)

    def execute(self, options, out_dir=None):
        """Executes a java JAR application."""
        # First we check if the jar file exists
        jar_file = self.get_jar_file()
        if not jar_file.endswith(".jar"):
            m = "{}: not a JAR file".format(jar_file)
            raise ProgramError(m)
        if not os.path.isfile(jar_file):
            m = "{}: no such file".format(jar_file)
            raise ProgramError(m)

        # Creating the required options
        self.required_options = self.get_child_required_options().copy()
        for key, value in self.get_jar_required_options().items():
            self.required_options[key] = value

        # Creating the command
        self.command = "{} {}".format(self.get_jar_command(),
                                      self.get_child_command())

        # Checking the options
        if "java_memory" not in options:
            options["java_memory"] = "4g"

        # Adding the jar file to the options
        options["jar_file"] = jar_file

        # Executing the jar file
        super(JAR, self).execute(options, out_dir)
