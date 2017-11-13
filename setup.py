from setuptools import setup, find_packages
import sys

import poreduck.version

if sys.version_info < (3,6):
    sys.exit('Sorry, Python < 3.6 is not supported')

with open("requirements.txt", 'r') as file:
    requirements = [line.strip() for line in file.readlines()]

long_description = """

Poreduck is a set of tools for handling, basecalling and producing quality metrics of MinION data.

"""

setup(
    name='poreduck',
    version=poreduck.__version__,
    packages=find_packages(),
    provides=['poreduck'],
    requires=['python (>=3.6)'],
    install_requires=requirements,
    url='github.com/alexiswl/poreduck',
    license='GPL',
    author='Alexis Lucattini',
    author_email='alexis.lucattini@agrf.org.au',
    description='Nanopore data handling poreduck.',
    package_dir={'poreduck': "poreduck"},
    entry_points={
        'console_scripts': ['poreduck=poreduck.poreduck_main:main']
    }
)

