import numpy as np
import pandas as pd
from tqdm import tqdm

from genome_tools import filter_df_to_interval, df_to_variant_intervals, VariantInterval

from genome_tools.data.extractors import TabixExtractor, VariantGenotypeExtractor
from genome_tools.data.get_variant_resolved_reads import extract_allelic_reads

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class VariantGenotypeLoader(PlotDataLoader):
    
    def _load(self, data: DataBundle, genotypes_vcf_path):
        with VariantGenotypeExtractor(genotypes_vcf_path) as extractor:
            variants = extractor[data.interval].rename(columns={'pos': 'end'})
        variants['start'] = variants['end'] - 1
        
        gt_mapping = {
            (0, 0): "A",
            (1, 1): "B",
        }
        
        variants["parsed_genotype"] = variants["gt"].map(gt_mapping)
        variants = variants.dropna(subset=["parsed_genotype"])

        data.variant_genotypes = variants
        return data


class GroupsByGenotypeLoader(PlotDataLoader):
    def _load(self, data: DataBundle, samples_metadata: pd.DataFrame, variant_interval: VariantInterval):
        ### Refactor to work with arbitrary groups
        variant_genotypes: pd.DataFrame = data.variant_genotypes # indiv_id, variant pairs

        # TODO: change filter_df_to_interval to accept VariantInterval
        variant_genotypes = filter_df_to_interval(variant_genotypes, variant_interval)
        if variant_genotypes.empty:
            raise ValueError("No genotypes found for the specified variant interval.")

        samples_with_genotype = samples_metadata.dropna(
            subset='indiv_id'
        ).reset_index(
            names='sample_id'
        ).merge(
            variant_genotypes
        ).set_index(
            "sample_id"
        ).rename(
            columns={'parsed_genotype': 'group'}
        ).sort_values(
            by="group"
        )
        data.groups_data = samples_with_genotype['group']
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


# TODO add fasta as input
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


class AllelicReadsLoaderFPTools(PlotDataLoader):

    def _load(self, data: DataBundle, sample_ids, samples_metadata: pd.DataFrame, variant_interval: VariantInterval):
        from footprint_tools.cutcounts import bamfile as BamFile
        if isinstance(sample_ids, (str, int, float)):
            sample_ids = [sample_ids]
        cram_paths = samples_metadata.loc[sample_ids, 'cram_file']
        reads = {}
        for sample_id, cram_path in tqdm(
            zip(sample_ids, cram_paths),
            total=len(sample_ids),
            desc="Allelic reads loader"
        ):
            extractor = BamFile(cram_path)
            # TODO replace with extractor
            reads[sample_id] = extractor.lookup_allelic(
                chrom=variant_interval.chrom,
                start=variant_interval.start,
                end=variant_interval.end,
                pos=variant_interval.pos,
                ref=variant_interval.ref,
                alt=variant_interval.alt
            )
        
        ref_cuts = np.zeros(len(data.interval))
        alt_cuts = np.zeros(len(data.interval))
        data.reads = reads
        for sample_id, allelic_reads in reads.items():

            sample_ref_cuts = allelic_reads[variant_interval.ref]["+"] + allelic_reads[variant_interval.ref]["-"]
            sample_alt_cuts = allelic_reads[variant_interval.alt]["+"] + allelic_reads[variant_interval.alt]["-"]
            # TODO: check the actual genotype
            if sum(sample_ref_cuts) == 0 or sum(sample_alt_cuts) == 0:
                print(f"Skipping sample {sample_id} due to no reads for one allele")
                continue
            ref_cuts += sample_ref_cuts
            alt_cuts += sample_alt_cuts
        data.ref_cuts = ref_cuts
        data.alt_cuts = alt_cuts
        data.variant_interval = variant_interval
        return data