#!/usr/bin/env python3

          
from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pybioutils',
    version='0.1',
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
    scripts=['src/cut_fasta.py', 'src/cut_list.py', 'src/new_python.py',
             ],
    packages=find_packages(exclude=["src/"]),
    include_package_data=True,
    install_requires=[],
    zip_safe=False,
)
