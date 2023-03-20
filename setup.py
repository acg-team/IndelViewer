########################################################################################################################
# Licensed Materials - Property of Gholamhossein Jowkar
#
#
# Copyright (C) 2023-2030 by Gholamhossein Jowkar
#
# This work was supported by the Swiss National Science Foundation (SNF) grants $31003A\_176316$ to Dr. M. A.
# The funding body did not play any role in the design of the study and collection, analysis, and interpretation of
# data and in writing the code.
#
########################################################################################################################
#
# This file is part of IndelViewer of ARPIP project: https://github.com/acg-team/IndelViewer
#
#                              ABOUT THE ARPIP PACKAGE
#                               =====================
# ARPIP: Ancestral Sequence Reconstruction with insertions and deletions under the Poisson Indel Process
# ARPIP is a joint maximum likelihood approach for phylogenetic ancestral sequence reconstruction, capable of modeling
# indels both biological and mathematically. You can find ARPIP in: https://github.com/acg-team/bpp-arpip
#
#
# IndelViewer is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# IndelViewer is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# You should have received a copy of the GNU Lesser General Public
# License along with IndeViewer. If not, see <http://www.gnu.org/licenses/>.
#
########################################################################################################################
# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: September 2022
# Modified by:
# Modified date:

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

with open('requirements.txt') as rf:
    requirements = rf.read()

setup(
    name='indelviewer',
    version='0.1.0',
    description='By using this library you can visualize indel events on the phylogeny tree.',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Gholam-Hossein Jowkar',
    author_email='jowk@zhaw.ch',
    url='https://github.com/acg-team/IndelViewer',
    license=license,
    packages=['indelviewer'],
    classifiers=[
            "Intended Audience :: Science/Research",
            "Intended Audience :: Developers",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.9",
            "Natural Language :: English",
            "Operating System :: MacOS",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX :: Linux",
            "Topic :: Software Development",
            "Topic :: Scientific/Engineering: Bioinformatics",
            "Topic :: Software Development :: Libraries :: Python Modules",
        ],
    install_requires=requirements
)