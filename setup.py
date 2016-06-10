from __future__ import absolute_import, division, print_function

import os
import sys

from setuptools import find_packages, setup

__version__ = "1.0"

install_requires = ["numpy>1.10", "scipy", "pysam>=0.8.2", "matplotlib"] 

setup(
	name = "genome_tools",
	version = __version__,
	description = "A toolkit for processing and plotting genomic datasets.",
	long_description = "",
	author = "Jeff Vierstra",
	author_email = "jvierstra@altiusinstitute.org",
	zip_safe = False,
	packages =  find_packages(),
    install_requires = install_requires,
)
