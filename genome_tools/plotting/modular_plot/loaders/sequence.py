from genome_tools.data.extractors import TabixExtractor, FastaExtractor
from genome_tools import df_to_genomic_intervals


from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class FastaLoader(PlotDataLoader):
    def _load(self, data: DataBundle, fasta_file):
        with FastaExtractor(fasta_file) as ext:
            data.sequence = ext[data.interval]
        return data


class MotifHitsLoader(PlotDataLoader):

    def _load(self, data: DataBundle, motif_annotations_path, motif_meta, best_by='dg'):
        with TabixExtractor(motif_annotations_path,
                             columns=['chrom', 'start', 'end', 'fp_id',
                                      'motif_chr', 'motif_start', 'motif_end',
                                      'pfm', 'dg', 'orient', 'sequence'
                             ]) as extractor:
            interval_motif_hits = extractor[data.interval]

        interval_motif_hits['dg'] = interval_motif_hits['dg'].astype(float)
        interval_motif_hits['start'] = interval_motif_hits['start'].astype(int)
        interval_motif_hits['end'] = interval_motif_hits['end'].astype(int)
        interval_motif_hits['motif_start'] = interval_motif_hits['motif_start'].astype(int)
        interval_motif_hits['motif_end'] = interval_motif_hits['motif_end'].astype(int)
        interval_motif_hits['motif_id'] = interval_motif_hits['pfm'].str.replace('.pfm', '')

        # FIXME: currently only keep the best hit per footprint
        if best_by == 'dg':
            interval_motif_hits = interval_motif_hits.groupby('fp_id', group_keys=False).filter(lambda x: x.nlargest(1, 'dg'))
        else:
            raise NotImplementedError(f'Unknown best_by method: {best_by}')

        interval_motif_hits = interval_motif_hits.merge(motif_meta, left_on='motif_id', right_index=True)
        data.motif_intervals = df_to_genomic_intervals(
            interval_motif_hits,
            data.interval,
            extra_columns=['orient', 'motif_start', 'motif_end', 'tf_name', 'pwm']
        )
        return data


