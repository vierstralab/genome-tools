from __future__ import absolute_import, division, print_function

import os
import sys

from setuptools import find_packages, setup

__version__ = "1.0.1"

install_requires = ["numpy>1.10", "scipy", "pysam>=0.8.2", "matplotlib"] 

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
    download_url = "https://github.com/jvierstra/genome-tools/archive/v1.0.1.tar.gz",
    classifiers=[
	    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
	    'Intended Audience :: Developers',      # Define that your audience are developers
	    'Topic :: Software Development :: Build Tools',
	    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
	    'Programming Language :: Python :: 2.7',     #Specify which pyhton versions that you want to support
	    'Programming Language :: Python :: 3',
	    'Programming Language :: Python :: 3.5', 
    ],
)
