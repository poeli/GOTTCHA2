#!/usr/bin/env python3

from setuptools import setup, find_packages
import os, subprocess, sys

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "gottcha2",
    version = "2.1.4",
    author = "Po-E Li",
    author_email = "poeli@lanl.gov",
    description = ("Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) v2 "),
    license = "BSD-3-Clause",
    keywords = ['bioinformatics', 'taxonomy', 'profiler'],
    url = "https://github.com/poeli/GOTTCHA2",
    packages=find_packages('GOTTCHA2'),
    install_requires=['numpy','pandas' ,'requests','tdqm','setuptools'],
    long_description=read('README.md'),
    entry_points={'console_scripts': ['gottcha2 = cmd:gottcha2_command',] },
    scripts=['gottcha2/gottcha2.py', 'gottcha2/pull_database.py', 'gottcha2/taxonomy.py', 'gottcha2/cmd.py']
)
