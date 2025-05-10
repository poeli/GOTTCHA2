#!/usr/bin/env python3

from setuptools import setup, find_packages
from gottcha.scripts.gottcha2 import __version__ as version

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name = "GOTTCHA2",
    version = version,
    author = "Po-E Li",
    author_email = "po-e@lanl.gov",
    description = ("Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) v2"),
    license = "BSD-3-Clause",
    keywords = ['bioinformatics', 'taxonomy', 'profiler', 'metagenomics'],
    url = "https://github.com/poeli/GOTTCHA2",
    packages=['gottcha'],
    package_dir={'GOTTCHA2': './'},
    install_requires=['numpy','pandas','requests','setuptools'],
    long_description=long_description,
    entry_points={'console_scripts': ['gottcha2 = cmd:gottcha2_command',] },
    scripts=['gottcha/scripts/gottcha2.py', 'gottcha/scripts/pull_database.py', 'gottcha/scripts/taxonomy.py', 'gottcha/scripts/cmd.py']
)
