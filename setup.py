#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip

import os
from setuptools import setup

from pgx_dna_seq import __version__

scripts = ["coverage_graph.py", "execute_pipeline.py", "automatic_report.py"]

setup(name="pgx_dna_seq",
      version=__version__,
      description="Automatic PGx DNA sequencing pipeline",
      author="Louis-Philippe Lemieux Perreault",
      author_email="louis-philippe.lemieux.perreault@statgen.org",
      url="http://www.statgen.org",
      license="GPL",
      scripts=[os.path.join("scripts", script_name) for script_name in scripts],
      install_requires=["numpy >= 1.8.1", "pandas >= 0.13.1", "ruffus >= 2.4.1", 
                        "matplotlib >=1.3.1"],
      packages=["pgx_dna_seq", "pgx_dna_seq.tool"],
      classifiers=['Operating System :: Linux',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3.4'])
