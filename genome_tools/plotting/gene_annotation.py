import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import numpy as np

import pysam

from genome_tools import genomic_interval

from .utils import pack_rows
from .colors.cm import get_vocab_color

ANNOTATION_KWDS = {
    'exon': {},
    'transcript': {},
    'utr': {}
}

def parse_gff_attribs(s):
    attr = {}
    for field in s.split(';'):
        (key, val) = field.split('=')
        attr[key] = val
    return attr

def load_gff_data(interval, gff_file):
    
    genes = {}
    transcripts = {}
    exons = {}
    cdss = {}
    utrs = {}

    fh = pysam.TabixFile(gff_file)
    
    for row in fh.fetch(interval.chrom, interval.start, interval.end, parser=pysam.asGTF()):
        
        feature_interval = genomic_interval(row.contig, row.start, row.end, strand = row.strand)
        attribs = parse_gff_attribs(row.attributes)

        if row.feature == 'gene':
            genes[attribs['gene_id']] = attribs['gene_name']
        elif row.feature == 'transcript':
           transcripts[attribs['gene_id']] = transcripts.get(attribs['gene_id'], []) + [(attribs['transcript_id'], feature_interval)]
        elif row.feature == 'exon':
            exons[attribs['transcript_id']] = exons.get(attribs['transcript_id'], []) + [(attribs['exon_id'], feature_interval)]
        elif row.feature == 'CDS':
            cdss[attribs['exon_id']] = cdss.get(attribs['exon_id'], []) + [(attribs['ID'], feature_interval)]
        elif row.feature == 'UTR':
            utrs[attribs['exon_id']] = utrs.get(attribs['exon_id'], []) + [(attribs['ID'], feature_interval)]

    fh.close()

    return genes, transcripts, exons, cdss, utrs


# Todo add color as styling kwargs
def gene_annotation_plot(interval, annotation_file, ax=None):

    if ax is None:
        ax = plt.gca()
    figure = ax.get_figure()

    genes, transcripts, exons, cdss, utrs = load_gff_data(interval, annotation_file)

    t = []
    for gene_id, gene_transcripts in transcripts.items():
        for transcript_id, transcript in gene_transcripts:
            t.append(transcript)

    row_indices = pack_rows(t)

    label_xoffset = len(interval) * 0.04

    label_arrow_height = 0.3
    label_arrow_width = len(interval) * 0.025
    label_arrow_xoffset = 0
    label_arrow_yoffset = 0.2

    for gene_id, gene_symbol in genes.items():
        
        if gene_id not in transcripts:
            continue

        for transcript_id, transcript in transcripts.get(gene_id, []):
            
            i = row_indices[transcript]

            over_left = over_right = start_visible = False

            if transcript.start < interval.start:
                over_left = True
                start = interval.start
            else:
                start = transcript.start

            if transcript.end > interval.end:
                over_right = True
                end = interval.end
            else:
                end = transcript.end

            # Draw trasnscript line
            ax.plot([start, end], [i, i], color='k', lw=2)

            for exon_id, exon in exons.get(transcript_id, []):
                
                for utr_id, utr in utrs.get(exon_id, []):
                    # make this a rectangle so we can put a border
                    l = mlines.Line2D([utr.start, utr.end], [i, i], color='k', lw=4)
                    ax.add_line(l)

                for cds_id, cds in cdss.get(exon_id, []):
                    l = mlines.Line2D([cds.start, cds.end], [i, i], color='blue', lw=5)
                    ax.add_line(l)

            # Draw labels/arrrow/etc...
            if transcript.strand == '+':

                label_arrow_xcoords = np.array([start, start, start + label_arrow_width])
                label_arrow_ycoords = np.array([i, i + label_arrow_height, i + label_arrow_height]) + label_arrow_yoffset
                ma = '>'

                if not over_left:
                    start_visible = True
                    label_xcoord = start + label_xoffset
                    label_ycoord = i + label_arrow_height + label_arrow_yoffset
                    ha = 'left'
                else:
                    label_xcoord = start - label_xoffset
                    label_ycoord = i
                    ha = 'right'

            else: 

                label_arrow_xcoords = np.array([end, end, end - label_arrow_width])
                label_arrow_ycoords = np.array([i, i + label_arrow_height, i + label_arrow_height]) + label_arrow_yoffset
                ma = '<'

                if not over_right:
                    start_visible = True
                    label_xcoord = end - label_xoffset
                    label_ycoord = i + label_arrow_height + label_arrow_yoffset
                    ha = 'right'
                else:
                    label_xcoord = end + label_xoffset
                    label_ycoord = i
                    ha = 'left'

            if start_visible:
                ax.plot(label_arrow_xcoords, label_arrow_ycoords, color = 'black')
                ax.plot(label_arrow_xcoords[-1], label_arrow_ycoords[-1], marker = ma, color = 'black', markeredgecolor='none')

            ax.text(label_xcoord, label_ycoord, gene_symbol, style = 'italic', verticalalignment = 'center', horizontalalignment = ha)

    ax.set_xlim(interval.start, interval.end)
    ax.set_ylim(bottom=-0.5, top=max(row_indices.values())+1)
    ax.yaxis.set_visible(False)