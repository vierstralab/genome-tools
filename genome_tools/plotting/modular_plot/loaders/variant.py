import numpy as np
import pandas as pd
from genome_tools.genomic_interval import filter_df_to_interval, df_to_variant_intervals, VariantInterval

from genome_tools.data.extractors import TabixExtractor, VariantGenotypeExtractor
from genome_tools.data.get_variant_resolved_reads import extract_allelic_reads

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class VariantGenotypeLoader(PlotDataLoader):
    
    def _load(self, data: DataBundle, vcf_path):
        with VariantGenotypeExtractor(vcf_path) as extractor:
            variants = extractor[data.interval].rename(columns={'pos': 'end'})
        variants['start'] = variants['end'] + 1
        
        gt_mapping = {
            (0, 0): "A",
            (1, 1): "B",
        }
        
        variants["parsed_genotype"] = variants["gt"].map(gt_mapping)
        variants = variants.dropna(subset=["parsed_genotype"])

        data.variant_genotype = variants
        return data


class FinemapLoader(PlotDataLoader):

    def _load(self, data: DataBundle, finemap_df: pd.DataFrame, region, trait, cs_id):
        finemap_df = finemap_df.query(
            f'region == "{region}" & trait == "{trait}" & cs_id == {cs_id}'
        ).drop_duplicates('end')
        data.unique_finemap_df = finemap_df
        return data


class AggregatedCAVLoader(PlotDataLoader):

    def _load(self, data: DataBundle, cavs_data, fdr_tr=0.1, color='k', notsignif_color='#C0C0C0'):
        filtered_cavs = filter_df_to_interval(cavs_data, data.interval)
        filtered_cavs['is_significant'] = filtered_cavs['min_fdr'] <= fdr_tr
        filtered_cavs['sig_es'] = np.clip(np.where(filtered_cavs['is_significant'], np.abs(filtered_cavs['logit_es_combined']), 0), 0, 2)
        group_ids_df = filtered_cavs.query('is_significant').groupby(['#chr', 'start', 'end', 'ref', 'alt'])['group_id'].apply(lambda x: ','.join(map(str, x))).reset_index()
        filtered_cavs = filtered_cavs.groupby(['#chr', 'start', 'end', 'ref', 'alt'], group_keys=False).apply(lambda x: x.nlargest(1, 'sig_es'))
        filtered_cavs = filtered_cavs.merge(group_ids_df, on=['#chr', 'start', 'end', 'ref', 'alt'], how='left', suffixes=('', '_list'))
        
        filtered_cavs['value'] = np.abs(filtered_cavs['logit_es_combined'])
        filtered_cavs['color'] = np.where(filtered_cavs['is_significant'], color, notsignif_color)
        data.cavs_intervals = df_to_variant_intervals(filtered_cavs, extra_columns=['value', 'color'])
        return data
    

class PerSampleCAVLoader(PlotDataLoader):
    __required_fields__ = ['nonaggregated_cavs_data']

    def _load(self, data: DataBundle, nonaggregated_cavs_data, sample_id, fdr_tr=0.1, color='k', notsignif_color='#C0C0C0'):
        filtered_cavs = TabixExtractor(nonaggregated_cavs_data)[data.interval].query(f'sample_id == "{sample_id}"')
        filtered_cavs['is_significant'] = filtered_cavs['FDR_sample'] <= fdr_tr
        filtered_cavs['sig_es'] = np.clip(np.where(filtered_cavs['is_significant'], np.abs(filtered_cavs['logit_es']), 0), 0, 2)
        
        filtered_cavs['value'] = np.abs(filtered_cavs['logit_es'])
        filtered_cavs['color'] = np.where(filtered_cavs['is_significant'], color, notsignif_color)
        data.cavs_intervals = df_to_variant_intervals(filtered_cavs, extra_columns=['value', 'color'])
        return data


class AllelicReadsLoader(PlotDataLoader):

    def _load(self, data: DataBundle, sample_ids, samples_metadata: pd.DataFrame, variant_interval: VariantInterval):
        if isinstance(sample_ids, (str, int, float)):
            sample_ids = [sample_ids]
        cram_paths = samples_metadata.loc[sample_ids, 'cram_file']
        reads = {}
        for sample_id, cram_path in zip(sample_ids, cram_paths):
            # TODO replace with extractor
            reads[sample_id] = extract_allelic_reads(cram_path, variant_interval, data.interval)
        data.reads = reads
        return data