# Copyright 2015 Jeff Vierstra

import numpy as np

class genomic_interval(object):
	"""Class that implements BED-style object
	"""

	def __init__(self, chrom, start, end, name = '.', score = None, strand = None):
		self.chrom = str(chrom)
		self.start = int(start)
		self.end = int(end)
		self.name = str(name)
		self.score = score
		self.strand = strand
		
	def __len__(self):
		"""returns length of element
		"""

		return self.end - self.start
	
	def __str__(self):
		"""returns a string-formated version of the element for printing
		"""
		return '\t'.join( [ str(x) for x in [self.chrom, self.start, self.end] ] )

	def widen(self, x, inplace=False):
		"""expands the coordinates

		:param x: number of nucletides to expand (can also be negative)
		:param inplace: modify estixing element or return a modified copy
		:type x: int
		:type inplace: bool
		:rtype: genomic_interval
		"""
		
		if inplace:
			self.start-=x
			self.end+=x
			return self
		return genomic_interval(self.chrom, self.start-x, self.end+x, self.name, self.strand)

	def shift(self, x, inplace=False):
		if inplace:
			self.start+=x
			self.end+=x
			return self
		else:	
			return genomic_interval(self.chrom, self.start+x, self.end+x, self.name, self.stand)

class genomic_interval_set(object):

	def __init__(self, iterator = []):
		self.intervals = []
		for interval in iterator:
			self.intervals.append(interval)
	
	def __len__(self):
		return len(self.intervals)

	def __iter__(self):
		for x in self.intervals:
			yield x

	def __getitem__(self, i):
		return self.intervals[i]

	def __iadd__(self, other):
		if type(other) == genomic_interval_set:
			self.intervals.extend(other.intervals)
		else:
			self.intervals.append(other)
		return self

def bin_intervals(intervals, bins):
	"""Bin genomic intervals by score thresholds

	Args:
		intervals (genomic_interval_set): a set of genomic intervals
						that contain a score value
		bins (array): a list of values comprising desired bins
						(must by monotonically increasing or decreasing)
	Returns:
		(list): a list of indices for each input interval to a corresponding bin
	"""
	scores = np.zeros(len(intervals))
	for i in range(len(intervals)):
		scores[i] = intervals[i].score
	return np.digitize(scores, thresh)


