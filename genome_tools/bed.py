# Copyright 2015 Jeff Vierstra

from . import GenomicInterval


def bed3_iterator(filehandle):
    """
    Generator that parses BED3 format from a string iterator

    Returns:
            genomic_interval
    """
    for line in filehandle:
        fields = line.strip().split()

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])

        yield GenomicInterval(chrom, start, end)


def bed5_iterator(filehandle):
    """
    Generator that parses BED5 format from a string iterator

    Returns:
            genomic_interval
    """

    for line in filehandle:
        fields = line.strip().split()

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        score = float(fields[4])

        yield GenomicInterval(chrom, start, end, name, score)


def bed6_iterator(filehandle):
    """
    Generator that parses BED6 format from a string iterator

    Returns:
            genomic_interval
    """

    for line in filehandle:
        fields = line.strip().split()

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        score = float(fields[4])
        strand = fields[5]

        yield GenomicInterval(chrom, start, end, name, score, strand)
