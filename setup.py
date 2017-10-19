#!/usr/bin/env python3

import sys          
from setuptools import setup, find_packages

__version__ = "0.1"


if sys.version_info < (3,5):
    print("python version must be at least 3.5")
    sys.exit(1)


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pybioutils',
    version='0.2',
    description='some things that I used',
    long_description=readme(),

    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        ],
    keywords='utilities bioinformatics',
    url='https://github.com/T-B-F/pyBioUtils/',
    author='Tristan Bitard-Feildel',
    author_email='tristan.bitard-feildel@impmc.upmc.fr',
    license='MIT',
    scripts=['BioUtils/bin/biotk',
             ],
    packages=find_packages(exclude=["BioUtils/bin/", "src/"]),
    include_package_data=True,
    install_requires=['biopython>=1.68',
                      'numpy',
                      'matplotlib'],

    zip_safe=False,
)
