"""pgx_dnaseq: A module to automatize and facilitate NGS data analysis."""

__author__ = "Louis-Philippe Lemieux Perreault and Abdellatif Daghrach"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault", "Abdellatif Daghrach",
               "Michal Blazejczyk"]
__license__ = "GPL"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


# Loading the verison
try:
    from .version import pgx_dnaseq_version as __version__
except ImportError:
    __version__ = None


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
