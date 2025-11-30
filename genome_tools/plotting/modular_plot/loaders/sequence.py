import pandas as pd
from typing import List
import numpy as np



from genome_tools import df_to_genomic_intervals, filter_df_to_interval, GenomicInterval, VariantInterval
from genome_tools.data.extractors import TabixExtractor, FastaExtractor
from genome_tools.data.pwm import read_pfm, get_allelic_scores, seq_logp, seq_heights_to_matrix

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class FastaLoader(PlotDataLoader):
    def _load(self, data: DataBundle, fasta_file):
        with FastaExtractor(fasta_file) as ext:
            data.sequence = ext[data.interval]
        return data
    

class OHESequenceLoader(PlotDataLoader):
    def _load(self, data: DataBundle):
        seq = data.sequence.upper()
        matrix = seq_heights_to_matrix(seq, np.ones(len(seq)))
        data.matrix = matrix
        return data


class MotifHitsLoader(PlotDataLoader):

    def _load(self, data: DataBundle, motif_annotations_path, motif_meta, min_motif_overlap=0.9):
        """
        Load motif hits overlapping annotation regions.
        
        Args:
            data (DataBundle): The data bundle containing annotation_regions.
            motif_annotations_path (str): Path to the tabix-indexed motif hits file (produced by MOODS) with columns ['#chr', 'start', 'end', 'motif_id', 'dg', 'orient', 'seq'].
            motif_meta (pd.DataFrame): DataFrame containing motif metadata indexed by motif_id.
            min_motif_overlap (float): Minimum fraction of motif overlap required to keep a hit.
        
        Returns:
            DataBundle: Updated data bundle with all_motifs_df.
        """

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
        return data


class MotifHitsSelectorLoader(PlotDataLoader):

    def _load(self, data, choose_by='dg', n_top_hits=1, motif_hits_threshold=None, variant_interval: VariantInterval=None):
        """
        Select motif hits based on a scoring metric and an optional selection threshold.

        Parameters
        ----------
        data : DataBundle
            Input data bundle containing `data.all_motifs_df` with all motif hits.

        choose_by : {'dg', 'ddg', 'concordant_ddg'}
            Scoring metric to use for selecting motif hits:
            - 'dg': Use the motif delta G score from the motif hits file.
            - 'weighted_dg': Weight the motif delta G score by per-position weights provided in data.seqeunce_weights
            - 'ddg': Compute absolute difference of motif scores between reference and alternate alleles.
            - 'concordant_ddg': Same as 'ddg', but keep only ddg values that match the direction of the variant effect (sign is taken from `variant_interval.value`).

        threshold : float or None
            Filter single best motif per region or filter by a cutoff:
            - None: Select the best motif per region
            - float: Keep all motifs with dg (or ddg) >= threshold

        variant_interval : VariantInterval, optional
            Required only for 'ddg' and 'concordant_ddg'. Concordance is determined based on the sign of `variant_interval.value`.

        Returns
        -------
        DataBundle
            Modified data bundle with `data.motif_intervals` containing the selected motif hits as genomic intervals.
        """
        motif_hits = data.all_motifs_df

        if choose_by == 'dg':
            metric_name = 'dg'
        
        elif choose_by == 'weighted_dg':
            metric_name = 'weighted_dg'
            
            if not hasattr(data, 'sequence_weights'):
                raise ValueError("data.sequence_weights is required for weighted_dg scoring. Attach another loader that sets sequence_weights to use this method.")
            assert len(data.sequence_weights) == len(data.interval), f"data.sequence_weights length must match data.interval length, got {len(data.sequence_weights)} vs {data.interval} ({len(data.interval)})"

            motif_hits = filter_df_to_interval(motif_hits, data.interval, strict=True)
            motif_hits['weighted_dg'] = motif_hits.apply(
                self.get_weighted_dg,
                axis=1,
                interval=data.interval,
                sequence_weights=data.sequence_weights,
            )

        elif choose_by in ('ddg', 'concordant_ddg'):
            if variant_interval is None:
                raise ValueError("variant_interval required for ddg scoring")
            if not isinstance(variant_interval, VariantInterval):
                raise ValueError("variant_interval must be an instance of genome_tools.VariantInterval")
            assert variant_interval.overlaps(data.interval), f"variant_interval must overlap data.interval, variant_interval={variant_interval.to_ucsc()} vs {data.interval}"
            metric_name = 'abs_ddg'

            motif_hits = filter_df_to_interval(motif_hits, variant_interval.gi())
            score_table = motif_hits.apply(
                self.score_row,
                axis=1,
                variant_interval=variant_interval
            )
            motif_hits[['ref_score', 'alt_score', 'pfm_matrix']] = score_table
            motif_hits['ddg'] = motif_hits['ref_score'] - motif_hits['alt_score']
            motif_hits['abs_ddg'] = motif_hits['ddg'].abs()

            if choose_by == 'concordant_ddg':
                effect_sign = np.sign(variant_interval.value)
                motif_hits = motif_hits.query(f"ddg * {effect_sign} > 0")

        else:
            raise ValueError(f"Unknown choose_by: {choose_by}")

        selected_hits = self._select_hits(motif_hits, metric_name, n_top_hits, motif_hits_threshold)

        data.motif_intervals = df_to_genomic_intervals(
            selected_hits,
            data.interval,
            extra_columns=['orient', 'region', 'tf_name', 'pfm_matrix']
        )
        return data
    
    def get_weighted_dg(self, row, interval: GenomicInterval, sequence_weights: np.ndarray):
        pfm_matrix = read_pfm(row['pfm'])
        weights = sequence_weights[
            row['start'] - interval.start: row['end'] - interval.end
        ]
        print(pfm_matrix.shape, len(weights))
        return seq_logp(
            mat=pfm_matrix,
            seq=row['seq'],
            weights=weights
        )

    @staticmethod
    def _select_hits(motif_hits: pd.DataFrame, metric_name: str, top: int, threshold: float | None):
        motif_hits = motif_hits.sort_values(metric_name, ascending=False)
        if threshold is not None:
            result = motif_hits.query(f"{metric_name} >= {threshold}")
        result = (
            motif_hits
            .groupby("region", group_keys=False)
            .head(top)
        )

        if 'pfm_matrix' not in result.columns:
            result['pfm_matrix'] = result['pfm'].apply(read_pfm)
        return result

    @staticmethod
    def score_row(row: pd.Series, variant_interval: VariantInterval) -> pd.Series:
        pfm_matrix = read_pfm(row['pfm'])
        offset = variant_interval.start - row['start']
        ref_score, alt_score = get_allelic_scores(
            pfm_matrix=pfm_matrix,
            sequence=row['seq'],
            ref=variant_interval.ref,
            alt=variant_interval.alt,
            offset=offset,
            orient=row['orient']
        )
        return pd.Series({
            "ref_score": ref_score,
            "alt_score": alt_score,
            "pfm_matrix": pfm_matrix
        })

