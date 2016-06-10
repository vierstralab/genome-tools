
# Copyright 2016 Jeff Vierstra

import numpy as np
import pysam

from . import genomic_interval

def load_tabix_data(file_path, interval, data_columns = [5]):
	""" Loads numeric data columns from a BED-format TABIX file """

	data = np.zeros((len(interval), len(data_columns)), dtype = np.float64)

	tabix = pysam.TabixFile(file_path)

	for row in tabix.fetch(interval.chrom, interval.start, interval.end, parser = pysam.asTuple()):
	    i = int(row[1])- interval.start
	    
	    for j, col in enumerate(data_columns):

	    	data[i, j] = row[col-1]

	tabix.close()
	    
	return data

