import re
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

import numpy as np

import pysam

from genome_tools.genomic_interval import GenomicInterval

from .utils import pack_rows
from .colors.cm import get_vocab_color

ANNOTATION_KWDS = {"exon": {}, "transcript": {}, "utr": {}}


def parse_gff_attribs(s):
    attr = {}
    for field in s.split(";"):
        if field.strip() == "":
            continue
        try:
            (key, val) = field.split("=")
        except ValueError:
            try:
                (key, val) = field.strip().split(" ")
            except:
                print(field, field.strip().split(" "))
                raise
        key = key.strip('"')
        val = val.strip('"')
        attr[key] = val
    return attr


def _delete_by_list(d, l):
    return [x[0] for k in l for x in d.pop(k, [])]


def load_gff_data(interval, gff_file, gene_symbol_exclude_regex=None):
    genes = {}
    transcripts = {}
    exons = {}
    cdss = {}
    utrs = {}

    fh = pysam.TabixFile(gff_file)

    for row in fh.fetch(
        interval.chrom, interval.start, interval.end, parser=pysam.asGTF()
    ):
        feature_interval = GenomicInterval(
            row.contig, row.start, row.end, strand=row.strand
        )
        attribs = parse_gff_attribs(row.attributes)

        if row.feature == "gene":
            genes[attribs["gene_id"]] = attribs["gene_name"]
        elif row.feature == "transcript":
            transcripts[attribs["gene_id"]] = transcripts.get(
                attribs["gene_id"], []
            ) + [(attribs["transcript_id"], feature_interval)]
        elif row.feature == "exon":
            exons[attribs["transcript_id"]] = exons.get(
                attribs["transcript_id"], []
            ) + [(attribs["exon_id"], feature_interval)]
        elif row.feature == "CDS":
            cdss[attribs["exon_id"]] = cdss.get(attribs["exon_id"], []) + [
                (attribs.get("ID"), feature_interval)
            ]
        elif row.feature == "UTR":
            utrs[attribs["exon_id"]] = utrs.get(attribs["exon_id"], []) + [
                (attribs.get("ID"), feature_interval)
            ]

    fh.close()

    gene_ids_to_delete = set()
    for gene_id, gene_symbol in genes.items():
        if gene_id not in transcripts:
            gene_ids_to_delete.add(gene_id)
        if gene_symbol_exclude_regex is not None and re.match(gene_symbol_exclude_regex, gene_symbol):
            gene_ids_to_delete.add(gene_id)
            transcript_ids = _delete_by_list(transcripts, genes[gene_id])
            exon_ids = _delete_by_list(exons, transcript_ids)
            _delete_by_list(cdss, exon_ids)
            _delete_by_list(utrs, exon_ids)

    for gene_id in gene_ids_to_delete:
        del genes[gene_id]

    return genes, transcripts, exons, cdss, utrs


# Todo add color as styling kwargs
def gene_annotation_plot(interval, annotation_file, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()

    genes, transcripts, exons, cdss, utrs = load_gff_data(interval, annotation_file, **kwargs)

    t = []
    for gene_id, gene_transcripts in transcripts.items():
        for transcript_id, transcript in gene_transcripts:
            t.append(transcript)

    row_indices = pack_rows(t)

    label_xoffset = len(interval) * 0.04

    label_arrow_height = 0.7
    label_arrow_width = len(interval) * 0.025
    label_arrow_xoffset = 0
    label_arrow_yoffset = 0.2

    for gene_id, gene_symbol in genes.items():

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
            ax.plot([start, end], [i, i], color="k", lw=0.75)

            for exon_id, exon in exons.get(transcript_id, []):
                for utr_id, utr in utrs.get(exon_id, []):
                    # make this a rectangle so we can put a border
                    l = mlines.Line2D([utr.start, utr.end], [i, i], color="k", lw=1.5)
                    ax.add_line(l)

                for cds_id, cds in cdss.get(exon_id, []):
                    l = mlines.Line2D([cds.start, cds.end], [i, i], color="blue", lw=2)
                    ax.add_line(l)

            # Draw labels/arrrow/etc...
            if transcript.strand == "+":
                label_arrow_xcoords = np.array(
                    [start, start, start + label_arrow_width]
                )
                label_arrow_ycoords = (
                    np.array([i, i + label_arrow_height, i + label_arrow_height])
                    + label_arrow_yoffset
                )
                ma = ">"

                if not over_left:
                    start_visible = True
                    label_xcoord = start + label_xoffset
                    label_ycoord = i + label_arrow_height + label_arrow_yoffset
                    ha = "left"
                else:
                    label_xcoord = start - label_xoffset
                    label_ycoord = i
                    ha = "right"

            else:
                label_arrow_xcoords = np.array([end, end, end - label_arrow_width])
                label_arrow_ycoords = (
                    np.array([i, i + label_arrow_height, i + label_arrow_height])
                    + label_arrow_yoffset
                )
                ma = "<"

                if not over_right:
                    start_visible = True
                    label_xcoord = end - label_xoffset
                    label_ycoord = i + label_arrow_height + label_arrow_yoffset
                    ha = "right"
                else:
                    label_xcoord = end + label_xoffset
                    label_ycoord = i
                    ha = "left"

            if start_visible:
                ax.plot(label_arrow_xcoords, label_arrow_ycoords, color="black", lw=0.5)
                ax.plot(
                    label_arrow_xcoords[-1],
                    label_arrow_ycoords[-1],
                    marker=ma,
                    color="black",
                    markeredgecolor="none",
                    ms=3,
                )

            ax.text(
                label_xcoord,
                label_ycoord,
                gene_symbol.replace('"', ''),
                style="italic",
                verticalalignment="center",
                horizontalalignment=ha,
                fontsize=5,
            )
            break

    ax.set_xlim(interval.start, interval.end)
    ax.set_ylim(bottom=-0.5, top=max(row_indices.values()) + 1)
    ax.yaxis.set_visible(False)
