import pandas as pd
from typing import List
import numpy as np

from genome_tools.data.extractors import TabixExtractor, FastaExtractor
from genome_tools import df_to_genomic_intervals, GenomicInterval


from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class FastaLoader(PlotDataLoader):
    def _load(self, data: DataBundle, fasta_file):
        with FastaExtractor(fasta_file) as ext:
            data.sequence = ext[data.interval]
        return data


class MotifHitsLoader(PlotDataLoader):

    def _load(self, data: DataBundle, motif_annotations_path, motif_meta, min_motif_overlap=0.9, best_by='dg'):

        annotation_regions: List[GenomicInterval] = [
            x for x in data.annotation_regions
        ]
        with TabixExtractor(
            motif_annotations_path,
            columns=['#chr', 'start', 'end', 'motif_id', 'dg', 'orient', 'seq']
        ) as extractor:
            regions_annotations = []
            for annotation_region in annotation_regions:
                interval_motif_hits = extractor[annotation_region]
                interval_motif_hits = interval_motif_hits
                interval_motif_hits['region'] = annotation_region.to_ucsc()
                interval_motif_hits['overlap'] = np.clip(
                    interval_motif_hits['end'],
                    annotation_region.start,
                    annotation_region.end
                ) - np.clip(
                    interval_motif_hits['start'],
                    annotation_region.start,
                    annotation_region.end
                )
                regions_annotations.append(interval_motif_hits)

        regions_annotations: pd.DataFrame = pd.concat(regions_annotations, ignore_index=True)
        regions_annotations['overlap_frac'] = regions_annotations.eval(
            'overlap / (end - start)'
        )
        # TMP FIX?
        regions_annotations['motif_id'] = regions_annotations['motif_id'].str.replace(
            '.pfm', ''
        )
        regions_annotations = regions_annotations.merge(
            motif_meta,
            left_on='motif_id',
            right_index=True
        )

        regions_annotations = regions_annotations.query(f'overlap_frac >= {min_motif_overlap}')

        data.all_motifs_df = regions_annotations

        filtered_regions = self.filter_motif_hits(
            regions_annotations,
            best_by=best_by
        )
        data.motif_intervals = df_to_genomic_intervals(
            filtered_regions,
            data.interval,
            extra_columns=['orient', 'region', 'tf_name', 'pfm']
        )
        return data

    def filter_motif_hits(self, motif_hits_df: pd.DataFrame, best_by='dg'):

        # FIXME: currently only keep the best hit per footprint
        if best_by == 'dg':
            motif_hits_df = (
                motif_hits_df
                .sort_values('dg', key=lambda s: s.abs(), ascending=False)
                .drop_duplicates(subset=['region'], keep='first')
            )
        else:
            raise NotImplementedError(f'Unknown best_by method: {best_by}')

        return motif_hits_df
