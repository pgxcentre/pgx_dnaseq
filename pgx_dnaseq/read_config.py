
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


import configparser

from . import tools
from .tools import *
from . import ProgramError


# Getting all the possible tools
_possible_tools = {}
for module_name in tools.__all__:
    for tool_name in eval(module_name).__all__:
        _possible_tools[tool_name] = eval("{}.{}".format(module_name,
                                                          tool_name))


def read_config_file(filename):
    """Reads a configuration file."""
    # Creating the configuration file parser (we want to read as string, so
    # that it is case sensitive)
    parser = configparser.RawConfigParser()
    parser.optionxform = str
    parser.read(filename)

    # Reading the configuration for each section and save it in a dict
    configuration = {}
    for section in parser.sections():
        options = {}
        for name, value in parser.items(section):
            options[name] = value
        configuration[section] = options

    # Returning the configuration
    return(configuration)


def get_pipeline_steps(filename):
    """Gets the different pipeline steps and their options."""
    # Will contain tuples (which tools to run with its their options)
    steps = []

    # Reading the configuration file
    pipeline = read_config_file(filename)

    # The keys should be steps number
    step_number = None
    try:
        step_number = sorted([int(i) for i in pipeline.keys()])
    except ValueError:
        m = "{}: {}: invalid steps".format(filename, list(pipeline.keys()))
        raise ProgramError(m)

    # Check that the step numbers are contiguous
    if step_number != list(range(step_number[0], step_number[-1] + 1)):
        m = "{}: {}: pipeline steps are not contiguous".format(filename,
                                                               step_number)
        raise ProgramError(m)

    # For each steps number, we gather the tool and its option
    for step in step_number:
        # Getting the step's tool and option
        tool_options = pipeline[str(step)]

        # Getting the tool name (and removing it from the options)
        if "tool" not in tool_options:
            m = "{}: step {}: no tool was specified".format(filename, step)
            raise ProgramError(m)
        tool_name = tool_options["tool"]
        if tool_name not in _possible_tools:
            m = "{}: {}: not a valid tool".format(filename, tool_name)
            raise ProgramError(m)
        del tool_options["tool"]

        # Saving the steps
        steps.append((_possible_tools[tool_name](), tool_options))

    # Returning the steps
    return steps
