# Created by:  Gholamhossein Jowkar <jowk@zhaw.ch>
# ACGT ZHAW
# Created date: September 2022
# Modified by:
# Modified date:

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='IndelViewer',
    version='0.1.0',
    description='Using this library you can visualize indel event on the phylogeny tree.',
    long_description=readme,
    author='Gholam-Hossein Jowkar',
    author_email='jowk@zhaw.ch',
    url='https://github.com/acg-team/arpip_indel_viewer',
    license=license,
    packages=['arpip.indelPoints.visualizer'],
    install_requires=['numpy','ete3','pandas']
)