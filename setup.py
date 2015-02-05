#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip

import os
from glob import glob
from setuptools import setup


MAJOR = 0
MINOR = 8
VERSION = "{}.{}".format(MAJOR, MINOR)


def write_version_file(fn=os.path.join("pgx_dnaseq", "version.py")):
    content = """
# THIS FILE WAS GENERATED AUTOMATICALLY BY PGX_DNASEQ SETUP.PY
pgx_dnaseq_version = {version}
"""
    a = open(fn, "w")
    try:
        a.write(content.format(version=VERSION))
    finally:
        a.close()


def setup_package():
    # Saving the version into a file
    write_version_file()

    setup(
        name="pgx_dnaseq",
        version=VERSION,
        description="Automatic PGx DNA sequencing pipeline",
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@statgen.org",
        url="http://www.statgen.org",
        license="CC BY-NC 4.0",
        scripts=glob(os.path.join("scripts", "*")),
        install_requires=["numpy >= 1.8.1", "pandas >= 0.13.1",
                          "ruffus >= 2.5", "matplotlib >=1.3.1",
                          "jinja2 >=2.7.3"],
        packages=["pgx_dnaseq", "pgx_dnaseq.tools"],
        package_data={"pgx_dnaseq": ["report_templates/*.tex",
                                     "report_templates/images/*"]},
        classifiers=['Operating System :: Linux',
                     'Programming Language :: Python',
                     'Programming Language :: Python :: 3.4'],
    )

    return


if __name__ == "__main__":
    setup_package()
