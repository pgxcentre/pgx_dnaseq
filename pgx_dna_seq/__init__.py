import configparser


class ProgramError(Exception):
    """An :py:class:`Exception` raised in case of a problem.
    
    :param msg: the message to print to the user before exiting.
    :type msg: string

    """
    def __init__(self, msg):
        """Construction of the :py:class:`ProgramError` class.

        :param msg: the message to print to the user
        :type msg: string

        """
        self.message = str(msg)

    def __str__(self):
        return self.message


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
