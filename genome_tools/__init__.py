"""
Toolkit for genomic intervals: read data for an interval (genome_tools.data.extractors),
manipulate intervals (genome_tools.genomic_interval) and plot them (genome_tools.plotting).

A guide for coding agents: AGENTS.md at the repository root
(https://github.com/vierstralab/genome-tools/blob/master/AGENTS.md).
"""
__all__ = ["bed", "GenomicInterval", "helpers"]

from .genomic_interval import *
