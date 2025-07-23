from __future__ import absolute_import, division, print_function

import os
import sys

from setuptools import find_packages, setup

__version__ = "2.0.0"

install_requires = ["numpy",
					"scipy",
					"pysam",
					"openpyxl",
					"pyBigWig",
					"matplotlib",
					"descartes",
					"shapely"] 

setup(
	name = "genome_tools",
	version = __version__,
	license = "GPL-3.0-or-later",
	description = "A toolkit for processing and plotting genomic datasets.",
	long_description = "",
	author = "Jeff Vierstra",
	author_email = "jvierstra@altius.org",
	zip_safe = False,
	packages =  find_packages(),
    install_requires = install_requires,
    download_url = "https://github.com/jvierstra/genome-tools/archive/v2.0.0.tar.gz",

)
