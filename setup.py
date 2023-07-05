#!/usr/bin/env python

from setuptools import setup, find_packages

with open('requirements.txt', 'r') as f:
      requirements = f.read().split('\n')

setup(name='pyreporter',
      version='0.2.2',
      description='Toolbox for finding reporter metabolites in genome-scale metabolic models using omics data',
      author='Roland Sauter',
      author_email='sauter.roland@gmail.com',
      url='',
      packages=find_packages('./pyreporter'),
      install_requires=requirements,
     )