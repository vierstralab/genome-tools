# Copyright 2015 Jeff Vierstra

import numpy as np

class genomic_interval(object):
	"""Class that implements BED-style object
	"""

	def __init__(self, chrom, start, end, name='.', score=None, strand=None, extra=None):
		self.chrom = str(chrom)
		self.start = int(start)
		self.end = int(end)
		self.name = str(name)
		self.score = score
		self.strand = strand
		self.extra = extra

	def __len__(self):
		"""Length of element
		"""
		return self.end - self.start
	
	def __str__(self):
		"""Returns a string-formated version of the element 
		for printing
		"""
		return '\t'.join( [ str(x) for x in [self.chrom, self.start, self.end] ] )

	def to_ucsc(self):
		"""Prints out a UCSC version of interval. Note
		that the start coordinate is shifted by 1 to
		reflect 1-based"""
		return f"{self.chrom}:{self.start}-{self.end}"

	def widen(self, x, inplace=False):
		"""Expands the coordinates
		"""
		if inplace:
			self.start -= x
			self.end += x
			return self
		return genomic_interval(self.chrom, self.start-x, self.end+x, self.name, self.strand)

	def shift(self, x, inplace=False):
		"""Shifts the coordinates"""
		if inplace:
			self.start += x
			self.end += x
			return self
		return genomic_interval(self.chrom, self.start+x, self.end+x, self.name, self.strand)

