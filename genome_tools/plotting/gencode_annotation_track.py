import numpy as np
import pysam

from matplotlib.patches import Rectangle

from .. import genomic_interval
from .track import track, pack_rows

class gencode_annotation_track(track):

	def __init__(self, interval,  **kwargs):
		super(gencode_annotation_track, self).__init__(interval, **kwargs)

		self.genes = {}
		self.transcripts = {}
		self.exons = {}
		self.cds = {}
		self.utrs = {}

	def parse_gff_attributes(self, s):
		attr = {}
		for field in s.split(';'):
			(key, val) = field.split('=')
			attr[key] = val
		return attr

	def load_data(self, filepath):
		"""Loads GFF data"""
		tabix = pysam.TabixFile(filepath)

		for row in tabix.fetch(self.interval.chrom, self.interval.start, self.interval.end, parser = pysam.asGTF()):

			feature_interval = genomic_interval(row.contig, row.start, row.end, strand = row.strand)
			a = self.parse_gff_attributes(row.attributes)

			if row.feature == "gene":
				self.genes[a["gene_id"]] = a["gene_name"]
			elif row.feature == "transcript":
				self.transcripts[a["gene_id"]] = self.transcripts.get(a["gene_id"], []) + [(a["transcript_id"], feature_interval)]
			elif row.feature == "exon":
				self.exons[a["transcript_id"]] = self.exons.get(a["transcript_id"], []) + [(a["exon_id"], feature_interval)]
			elif row.feature == "CDS":
				self.cds[a["exon_id"]] = self.cds.get(a["exon_id"], []) + [(a["ID"], feature_interval)]
			elif row.feature == "UTR":
				self.utrs[a["exon_id"]] = self.utrs.get(a["exon_id"], []) + [(a["ID"], feature_interval)]

		tabix.close()
	
	def get_transcript_intervals(self):
		transcript_intervals = []
		for gene, trans in self.transcripts.items():
			for transcript_name, transcript_interval in trans:
				transcript_intervals.append(transcript_interval)

		return transcript_intervals

	def format_axis(self, ax):

		super(gencode_annotation_track, self).format_axis(ax)

		ax.xaxis.set(visible = False)
		ax.yaxis.set(visible = False)
		
	def render(self, ax, padding = 5):

		self.format_axis(ax)
		self.format_spines(ax, remove_spines=['top', 'bottom', 'left', 'right'])
		# pack transcripts

		trans = self.get_transcript_intervals()
		nrows, rows = pack_rows(trans, padding = padding)

		label_xoffset = len(self.interval) * 0.04

		label_arrow_height = 0.3
		label_arrow_width = len(self.interval) * 0.025

		label_arrow_xoffset = 0
		label_arrow_yoffset = 0.2

		for gene, gene_name in self.genes.items():
	
			if not gene in self.transcripts:
				continue
	
			for (transcript, transcript_interval) in self.transcripts[gene]:
			
				n = rows[transcript_interval]

				x0 = transcript_interval.start #- self.interval.start
				x1 = transcript_interval.end  #- transcript_interval.start + x0
				y = n

				over_left = x0 < self.interval.start #0
				over_right = x1 > self.interval.end #len(self.interval)
				start_visible = False
			
				ax.plot([max(x0, self.interval.start), min(x1, self.interval.end)], [y, y], color = 'grey', lw = 2, zorder = 0)

				if transcript_interval.strand == '+':

					label_arrow_xcoords = np.array([x0, x0, x0 + label_arrow_width])
					label_arrow_ycoords = np.array([y, y + label_arrow_height, y + label_arrow_height ]) + label_arrow_yoffset
					ma = '>'

					if not over_left:

						start_visible = True

						label_xcoord = x0 + label_xoffset
						label_ycoord = y + label_arrow_height + label_arrow_yoffset
						ha = "left"

					else:

						label_xcoord = self.interval.start - label_xoffset
						label_ycoord = y
						ha = "right"

				else: 

					label_arrow_xcoords = np.array([x1, x1, x1 - label_arrow_width])
					label_arrow_ycoords = np.array([y, y + label_arrow_height, y + label_arrow_height ]) + label_arrow_yoffset
					ma = '<'

					if not over_right:

						start_visible = True

						label_xcoord = x1 - label_xoffset
						label_ycoord = y + label_arrow_height + label_arrow_yoffset
						ha = "right"

					else:

						label_xcoord = self.interval.end + label_xoffset
						label_ycoord = y
						ha = "left"

				if start_visible:
					ax.plot(label_arrow_xcoords, label_arrow_ycoords, color = 'black')
					ax.plot(label_arrow_xcoords[-1], label_arrow_ycoords[-1], marker = ma, color = 'black', markeredgecolor='none')

				ax.text(label_xcoord, label_ycoord, gene_name, style = 'italic', verticalalignment = 'center', horizontalalignment = ha)
				
				if over_left:
					ax.plot(self.interval.start, y, marker = ma, clip_on = False, color = 'grey', markeredgecolor='none')
				if over_right:
					ax.plot(self.interval.end, y, marker = ma, clip_on = False, color = 'grey', markeredgecolor='none')
					
				# Loop through exons
				for (exon, exon_interval) in self.exons.get(transcript, []):

					for (utr, utr_interval) in self.utrs.get(exon, []):
						p = Rectangle((utr_interval.start, -0.1+n), utr_interval.end-utr_interval.start, 0.2, edgecolor = 'none', facecolor = 'goldenrod', zorder = 1)
						ax.add_patch(p)

					for (cds, cds_interval) in self.cds.get(exon, []):
						p = Rectangle((cds_interval.start, -0.15+n), cds_interval.end-cds_interval.start, 0.3, edgecolor = 'none', facecolor = 'goldenrod', zorder = 2)
						ax.add_patch(p)
	
		ax.set_ylim(bottom = 0.5, top = nrows + 1.5)
		
		super(gencode_annotation_track, self).render(ax)
		
