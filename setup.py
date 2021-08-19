#!/usr/bin/env python3






from distutils.core import setup

setup(name='gottcha',
      version='2.1.4',
      description='Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)',
      author='Po-E Li',
      author_email='poeli@lanl.gov',
      url='https://github.com/poeli/GOTTCHA2',
      packages=['gottcha2'],
      package_dir={'gottcha2': './'},
      scripts=['gottcha2.py', 'pull_database.py', 'taxonomy.py']
     )
