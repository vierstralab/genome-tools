import pandas as pd
from typing import List

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

        motif_regions: List[GenomicInterval] = [
            x for x in data.annotation_regions
        ]
        with TabixExtractor(
            motif_annotations_path,
            columns=[]
        ) as extractor:
            regions_annotations = []
            for motif_region in motif_regions:
                interval_motif_hits = extractor[motif_region]
                interval_motif_hits = interval_motif_hits.rename(
                    columns={
                        'start': 'motif_start',
                        'end': 'motif_end'
                    }
                )
                interval_motif_hits['start'] = motif_region.start
                interval_motif_hits['end'] = motif_region.end
                interval_motif_hits['overlap_frac'] = (
                    interval_motif_hits[['motif_end', 'end']].min(axis=1) -
                    interval_motif_hits[['motif_start', 'start']].max(axis=1)
                ) + 1 / len(motif_region)
                regions_annotations.append(interval_motif_hits)

        # TMP FIX
        interval_motif_hits['motif_id'] = interval_motif_hits['pfm'].str.replace(
            '.pfm', ''
        )
        interval_motif_hits = interval_motif_hits.merge(
            motif_meta,
            left_on='motif_id',
            right_index=True
        )

        interval_motif_hits = self.filter_motif_hits(
            interval_motif_hits,
            min_motif_overlap=min_motif_overlap,
            best_by=best_by
        )
        data.motif_intervals = df_to_genomic_intervals(
            interval_motif_hits,
            data.interval,
            extra_columns=['orient', 'motif_start', 'motif_end', 'tf_name', 'pwm']
        )
        return data

    def filter_motif_hits(self, motif_hits_df: pd.DataFrame, min_motif_overlap=0.9, best_by='dg'):
        motif_hits_df = motif_hits_df.query(f'overlap_frac >= {min_motif_overlap}')

        # FIXME: currently only keep the best hit per footprint
        if best_by == 'dg':
            motif_hits_df = motif_hits_df.groupby(
                ['region_start', 'region_end'],
                group_keys=False
            ).filter(lambda x: x.nlargest(1, 'dg'))
        else:
            raise NotImplementedError(f'Unknown best_by method: {best_by}')

        return motif_hits_df
