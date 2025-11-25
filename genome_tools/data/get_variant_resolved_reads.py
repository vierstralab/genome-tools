import pysam

import numpy as np

from .extractors import OverlappingReadsExtractor
from genome_tools import VariantInterval, GenomicInterval


class GenotypeError(Exception):
	pass

class DiploidError(Exception):
	pass

class ReadBiasError(Exception):
	pass

class ReadAlignmentError(Exception):
	pass

class ReadGenotypeError(Exception):
	pass


def get_base_call(pileupread):
	if pileupread is None:
		return None
	if pileupread.query_position is None:
		return None
	else:
		return pileupread.alignment.query_sequence[pileupread.query_position]


def check_alleles(pileupread, ref_allele, nonref_allele):

	if pileupread is None:
		return True

	# if pileupread.alignment.mapping_quality<30:
	# 	raise ReadAlignmentError()
	
	read_allele = get_base_call(pileupread)
	if read_allele != ref_allele and read_allele != nonref_allele:
		return ReadGenotypeError()


def get_5p_offset(pileupread):
	"""
	Returns position of variant relative to 5' of read
	"""
	if pileupread.query_position is None: # pileup overlaps deletion 
		return None
	elif pileupread.alignment.is_reverse:
		return pileupread.alignment.query_length-pileupread.query_position
	else:
		return pileupread.query_position+1


def check_bias(pileupread, offset=3, baseq=20):

	if pileupread is None:
		return True

	if get_5p_offset(pileupread)<=offset:
		raise ReadBiasError()

	# if get_base_quality(pileupread)<baseq:
	# 	raise ReadBiasError()

	return True


def check_read(read, ref, alt):
    check_alleles(read, ref, alt) # only returns true if read exists
    check_bias(read) # only returns true if read exists
    read_allele = get_base_call(read) # returns None if read doesn't exist
    return read_allele

def check_reads(reads_1, reads_2, unique_reads, ref, alt):
    for read in unique_reads:
        try:
            read1 = reads_1.get(read, None)
            read1_allele = check_read(read1, ref, alt)

            read2 = reads_2.get(read, None)
            read2_allele = check_read(read2, ref, alt)
            
            # read is either first or second in the pair
            read_allele = read1_allele or read2_allele

        except (ReadBiasError, ReadGenotypeError):
            continue

        yield read_allele, (read1 or read2).alignment


def choose_reads_for_plottng(assigned_reads, reads1, reads2, variant_interval, seed=None):
	"""
	Choose reads for plotting
	"""
	reads = []
	rng = np.random.default_rng(seed=seed)
	for read_name, read_base in assigned_reads:
		if read_name not in reads1:
			read = reads2[read_name].alignment
		elif read_name not in reads2:
			read = reads1[read_name].alignment
		else: # both mates are overlapping the interval
			first = reads1[read_name].alignment
			second = reads2[read_name].alignment
			if first.reference_start <= variant_interval.start and first.reference_end >= variant_interval.end:
				read = first
			elif second.reference_start <= variant_interval.start and second.reference_end >= variant_interval.end:
				read = second
			else:
				read = rng.choice([first, second])
		reads.append(
               GenomicInterval(read.reference_name, read.reference_start, read.reference_end, base=read_base)
        )

	return reads


def extract_allelic_reads(
        bam_file_path,
        variant_interval: VariantInterval,
        genomic_interval: GenomicInterval=None
	):
	if genomic_interval is None:
		genomic_interval = variant_interval
	assert genomic_interval.chrom == variant_interval.chrom
	assert genomic_interval.start <= variant_interval.start
	assert genomic_interval.end >= variant_interval.end
	
	with OverlappingReadsExtractor(bam_file_path) as extractor:
		reads_1, reads_2, read_pairs = extractor[genomic_interval]
		variant_reads_1, variant_reads_2, variant_read_pairs = extractor[variant_interval]
	
	assigned_reads = [(x, "N") for x in read_pairs if x not in variant_read_pairs]
	for read_base, read in check_reads(
		variant_reads_1, variant_reads_2, variant_read_pairs, variant_interval.ref, variant_interval.alt):
		assigned_reads.append((read.query_name, read_base))

	return choose_reads_for_plottng(assigned_reads, reads_1, reads_2, variant_interval)
