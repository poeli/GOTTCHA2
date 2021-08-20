#!/usr/bin/env python3





from setuptools import setup

setup(name='gottcha',
      version='2.1.4',
      description='Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)',
      author='Po-E Li',
      author_email='poeli@lanl.gov',
      url='https://github.com/poeli/GOTTCHA2',
      packages=['gottcha2'],
      package_dir={'gottcha2': './'},
      scripts=['gottcha2.py', 'pull_database.py', 'taxonomy.py', 'cmd.py'],
      install_requires=['minimap2>=2.1','gawk','pandas','biom-format>=2.1.7','requests','tdqm'],
      entry_points={ 'console_scripts': [ 'gottcha2 = gottcha2.cmd:gottcha2_command']}
)
