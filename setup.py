#!/usr/bin/env python

from setuptools import setup

setup(
    name='fampipeline',
    version='0.0.1',
    author='Bernie Pope, Khalid Mahmood, Jessica Chung',
    author_email='kmahmood@unimelb.edu.au',
    packages=['src'],
    entry_points={
        'console_scripts': ['fampipeline = src.main:main']
    },
    url='https://github.com/khalidm/fampipeline',
    license='LICENSE.txt',
    description='fampipeline is a pipeline system for bioinformatics workflows\
     with support for running pipeline stages on a distributed compute cluster.',
    long_description=open('README.md').read(),
    install_requires=[
        "ruffus",
        "drmaa",
        "PyYAML"
    ],
)
