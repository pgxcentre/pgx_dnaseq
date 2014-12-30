#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip

import os
from glob import glob
from setuptools import setup

from pgx_dna_seq import __version__

setup(name="pgx_dna_seq",
      version=__version__,
      description="Automatic PGx DNA sequencing pipeline",
      author="Louis-Philippe Lemieux Perreault",
      author_email="louis-philippe.lemieux.perreault@statgen.org",
      url="http://www.statgen.org",
      license="GPL",
      scripts=glob(os.path.join("scripts", "*")),
      install_requires=["numpy >= 1.8.1", "pandas >= 0.13.1", "ruffus >= 2.4.1", 
                        "matplotlib >=1.3.1", "jinja2 >=2.7.3"],
      packages=["pgx_dna_seq", "pgx_dna_seq.tool"],
      package_data={"pgx_dna_seq": ["report_templates/*.tex",
                                    "report_templates/images/*"]},
      classifiers=['Operating System :: Linux',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3.4'])
