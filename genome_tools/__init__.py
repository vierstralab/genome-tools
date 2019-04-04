__all__ = ["bed", "genomic_interval", "helpers"]

from .genomic_interval import *

import pysam

def load_data(filepath, interval, data_columns=[5], dtype=np.float):

	""" Loads numeric data columns from a BED-format TABIX file """

	res = np.zeros((len(interval), len(data_columns)), dtype = dtype)

	tabix = pysam.TabixFile(filepath)

	for row in tabix.fetch(interval.chrom, interval.start, interval.end, parser = pysam.asTuple()):
	    i = int(row[1])- interval.start
	    
	    for j, col in enumerate(data_columns):
	    	res[i, j] = dtype(row[col-1])

	tabix.close()
	    
	return res


def load_segmented_data(filepath, interval):
	
	res = genomic_interval_set()

	tabix = pysam.TabixFile(filepath)
	for row in tabix.fetch(interval.chrom, interval.start, interval.end, parser = pysam.asTuple()):
		
		chrom=row[0]
		start = int(row[1])
		end = int(row[2])
		
		try:
			name=row[3]
		except:
			name='.'

		try:
			score=float(row[4])
		except:
			score=-np.inf
	
		try:
			strand=row[5]
		except:
			strand='+'

		res += genomic_interval(chrom, start, end, name=name, score=score, strand=strand) 
	
	tabix.close()

	return res

