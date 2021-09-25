#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='pyreporter',
      version='0.2',
      description='Algorithm for finding reporter metabolites in a GEM using omics data',
      author='Roland Sauter',
      author_email='sauter.roland@gmail.com',
      url='',
      packages=find_packages('./pyreporter'),
      install_requires=[
            'numpy', 'pandas', 'cobra',
            'scikit-network', 'networkx',
            ],
     )