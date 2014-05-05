#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip

import os
from setuptools import setup

setup(name="pgx_dna_seq",
      version="0.1",
      description="Automatic PGx DNA sequencing pipeline",
      author="Louis-Philippe Lemieux Perreault",
      author_email="louis-philippe.lemieux.perreault@statgen.org",
      url="http://www.statgen.org",
      license="GPL",
##       scripts=[os.path.join("scripts", "manhattan_generator"),],
##       install_requires=["matplotlib >=1.3.1", "numpy >= 1.8.0"],
      classifiers=['Operating System :: Linux',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3.4'])
