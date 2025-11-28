import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import List

from genome_tools import filter_df_to_interval, df_to_variant_intervals, VariantInterval, GenomicInterval

from genome_tools.data.extractors import TabixExtractor, VariantGenotypeExtractor
from genome_tools.data.get_variant_resolved_reads import extract_allelic_reads

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class VariantGenotypeLoader(PlotDataLoader):
    
    def _load(self, data: DataBundle, genotypes_vcf_path, samples_metadata: pd.DataFrame):
        with VariantGenotypeExtractor(genotypes_vcf_path) as extractor:
            variants = extractor[data.interval].rename(columns={'pos': 'end'})
        variants['start'] = variants['end'] - 1
        
        gt_mapping = {
            (0, 0): "AA",
            (1, 1): "BB",
            (0, 1): "AB",
            (1, 0): "AB",
        }
        
        variants["parsed_genotype"] = variants["gt"].map(gt_mapping)
        variants = variants.dropna(subset=["parsed_genotype"])

        variants_with_sample_id = samples_metadata.dropna(
            subset='indiv_id'
        ).reset_index(
            names='sample_id'
        ).merge(
            variants
        ).set_index(
            "sample_id"
        )

        data.variant_genotypes = variants_with_sample_id
        return data


class GroupsByGenotypeLoader(PlotDataLoader):
    def _load(self, data: DataBundle, variant_interval: VariantInterval, groups=["AA", "BB"]):
        ### Refactor to work with arbitrary groups
        variant_genotypes: pd.DataFrame = data.variant_genotypes # indiv_id, variant pairs

        # TODO: change filter_df_to_interval to accept VariantInterval
        variant_genotypes = filter_df_to_interval(variant_genotypes, variant_interval)
        if variant_genotypes.empty:
            raise ValueError("No genotypes found for the specified variant interval.")
        variant_genotypes = variant_genotypes.rename(
            columns={'parsed_genotype': 'group'}
        )
        variant_genotypes = variant_genotypes[
            variant_genotypes['group'].isin(groups)
        ].sort_values(
            by="group"
        )

        data.groups_data = variant_genotypes['group']
        return data

class FinemapLoader(PlotDataLoader):

    def _load(self, data: DataBundle, finemap_df: pd.DataFrame, region, trait, cs_id):
        finemap_df = finemap_df.query(
            f'region == "{region}" & trait == "{trait}" & cs_id == {cs_id}'
        ).drop_duplicates('end')
        data.unique_finemap_df = finemap_df
        return data


class AggregatedCAVLoader(PlotDataLoader):

    def _load(self, data: DataBundle, aggregated_cavs_df, fdr_tr=0.1, color='k', choose_best_by='sig_es', notsignif_color='#C0C0C0'):
        filtered_cavs = filter_df_to_interval(aggregated_cavs_df, data.interval)
        filtered_cavs['is_significant'] = filtered_cavs['min_fdr'] <= fdr_tr
        filtered_cavs['sig_es'] = np.clip(np.where(filtered_cavs['is_significant'], np.abs(filtered_cavs['logit_es_combined']), 0), 0, 2)
        significant_groups = filtered_cavs.query('is_significant').groupby(
            ['#chr', 'start', 'end', 'ref', 'alt']
        )['group_id'].apply(
            lambda x: ','.join(map(str, x))
        )
        if choose_best_by == "sig_es":
            filtered_cavs = (
                filtered_cavs
                .sort_values('sig_es', key=lambda s: s.abs(), ascending=False)
                .drop_duplicates(subset=['#chr', 'start', 'end', 'ref', 'alt'], keep='first')
            )
        filtered_cavs = filtered_cavs.set_index(['#chr', 'start', 'end', 'ref', 'alt'])
        filtered_cavs['significant_groups'] = significant_groups
        
        filtered_cavs['value'] = np.abs(filtered_cavs['logit_es_combined'])
        filtered_cavs['color'] = np.where(filtered_cavs['is_significant'], color, notsignif_color)
        data.cavs_intervals = df_to_variant_intervals(
            filtered_cavs.reset_index(),
            extra_columns=['value', 'is_significant', 'color', 'significant_groups']
        )
        return data
    

class PerSampleCAVLoader(PlotDataLoader):
    __required_fields__ = ['nonaggregated_cavs_data']

    def _load(self, data: DataBundle, nonaggregated_cavs_tabix, sample_id, fdr_tr=0.1, color='k', notsignif_color='#C0C0C0'):
        with TabixExtractor(nonaggregated_cavs_tabix) as extractor:
            filtered_cavs = extractor[data.interval].query(f'sample_id == "{sample_id}"')

        filtered_cavs['is_significant'] = filtered_cavs['FDR_sample'] <= fdr_tr
        filtered_cavs['sig_es'] = np.clip(np.where(filtered_cavs['is_significant'], np.abs(filtered_cavs['logit_es']), 0), 0, 2)
        
        filtered_cavs['value'] = np.abs(filtered_cavs['logit_es'])
        filtered_cavs['color'] = np.where(filtered_cavs['is_significant'], color, notsignif_color)
        data.cavs_intervals = df_to_variant_intervals(
            filtered_cavs,
            extra_columns=['value', 'color']
        )
        return data


class AllelicReadsLoader(PlotDataLoader):

    def _load(self, data: DataBundle, samples_metadata: pd.DataFrame, variant_interval: VariantInterval, sample_ids=None):
        assert variant_interval.overlaps(data.interval), f"variant_interval must overlap data.interval. Got {variant_interval.to_str()} and {data.interval.to_ucsc()}"

        variant_data: pd.DataFrame = data.variant_genotypes.query('parsed_genotype == "AB"')
        variant_data = filter_df_to_interval(variant_data, variant_interval)
        if variant_data.empty:
            raise ValueError("No heterozygous genotypes found for the specified variant_interval.")
        assert variant_data.index.is_unique

        if sample_ids is None:
            sample_ids = variant_data.index.tolist()
        elif isinstance(sample_ids, (str, int, float)):
            sample_ids = [sample_ids]

        cram_paths = samples_metadata.loc[sample_ids, 'cram_file']

        reads = {}
        for sample_id, cram_path in tqdm(zip(sample_ids, cram_paths), total=len(sample_ids), desc="Allelic reads loader"):
            extracted_reads = extract_allelic_reads(cram_path, variant_interval, data.interval)
            if not self.check_reads(extracted_reads, variant_interval):
                continue
            reads[sample_id] = extracted_reads
        data.reads = reads
        return data
    
    def check_reads(self, extracted_reads: List[GenomicInterval], variant_interval: VariantInterval):
        n_ref, n_alt = 0, 0
        for read in extracted_reads:
            n_ref += read.base == variant_interval.ref
            n_alt += read.base == variant_interval.alt
        if n_ref > 0 and n_alt > 0:
            return True
        return False


class AllelicReadsLoaderFPTools(PlotDataLoader):

    def _load(self, data: DataBundle, samples_metadata: pd.DataFrame, variant_interval: VariantInterval, sample_ids=None):
        from footprint_tools.cutcounts import BamFileExtractor
        assert variant_interval.overlaps(data.interval), f"variant_interval must overlap data.interval. Got {variant_interval.to_str()} and {data.interval.to_ucsc()}"

        variant_data: pd.DataFrame = data.variant_genotypes.query('parsed_genotype == "AB"')
        variant_data = filter_df_to_interval(variant_data, variant_interval)
        if variant_data.empty:
            raise ValueError("No heterozygous genotypes found for the specified variant_interval.")
        assert variant_data.index.is_unique

        if sample_ids is None:
            sample_ids = variant_data.index.tolist()
        elif isinstance(sample_ids, (str, int, float)):
            sample_ids = [sample_ids]

        cram_paths = samples_metadata.loc[sample_ids, 'cram_file']
        reads = {}
        for sample_id, cram_path in tqdm(
            zip(sample_ids, cram_paths),
            total=len(sample_ids),
            desc="Allelic reads loader fp-tools"
        ):
            extractor = BamFileExtractor(cram_path)
            reads[sample_id] = extractor.lookup_allelic(
                chrom=variant_interval.chrom,
                start=data.interval.start,
                end=data.interval.end,
                pos=variant_interval.pos - 1, # weird 0-based in fp-tools
                ref=variant_interval.ref,
                alt=variant_interval.alt
            )
            extractor.close()
        
        data.ref_reads = self.convert_reads_to_list(reads, variant_interval.ref)
        data.alt_reads = self.convert_reads_to_list(reads, variant_interval.alt)
        data.variant_interval = variant_interval
        return data
    
    @staticmethod
    def convert_reads_to_list(reads_dict, allele: str):
        reads = []
        for sample_id, sample_reads_dict in reads_dict.items():
            for read in sample_reads_dict[allele]['fragments']:
                assert hasattr(read, 'is_reverse')
                read.sample_id = sample_id
                read.base = allele
                reads.append(read)
        return reads