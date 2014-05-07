import configparser

from pgx_dna_seq import tool
from pgx_dna_seq.tool import *
from pgx_dna_seq import ProgramError


# Getting all the possible tools
__possible_tools = {}
for module_name in tool.__all__:
    for tool_name in eval(module_name).__all__:
        __possible_tools[tool_name] = eval("{}.{}".format(module_name,
                                                          tool_name))

def read_config_file(filename):
    """Reads a configuration file."""
    # Creating the configuration file parser (we want to read as string, so that
    # it is case sensitive)
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
        if tool_name not in __possible_tools:
            m = "{}: {}: not a valid tool".format(filename, tool_name)
            raise ProgramError(m)
        del tool_options["tool"]

        # Saving the steps
        steps.append((__possible_tools[tool_name](), tool_options))

    # Returning the steps
    return steps