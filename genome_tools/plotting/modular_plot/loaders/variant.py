import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import List

from genome_tools import filter_df_to_interval, df_to_variant_intervals, VariantInterval, GenomicInterval

from genome_tools.data.extractors import TabixExtractor, VariantGenotypeExtractor, AllelicReadsExtractor

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class VariantIntervalLoader(PlotDataLoader):
    required_loader_kwargs = ['variant_interval']


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
    def _load(self, data: DataBundle, groups=["AA", "BB"]):
        variant_interval: VariantInterval = data.variant_interval
        interval: GenomicInterval = data.interval
        variant_genotypes: pd.DataFrame = data.variant_genotypes # sample_id, variant pairs

        assert variant_interval.overlaps(interval), f"variant_interval must overlap data.interval. Got {variant_interval.to_str()} and {interval.to_ucsc()}"
        # groups
        variant_genotypes = filter_df_to_interval(variant_genotypes, variant_interval)
        if variant_genotypes.empty:
            raise ValueError("No genotypes found for the specified variant interval.")

        variant_genotypes['group'] = variant_genotypes['parsed_genotype'].map(
            {
                g: f"group{i}" for i, g in enumerate(groups, 1)
            }
        )
        filtered_genotypes = variant_genotypes.dropna(subset=['group'])

        if filtered_genotypes.empty:
            gts = variant_genotypes['parsed_genotype'].unique().tolist()
            raise ValueError(f"No samples with GT={groups} in interval. Available genotypes: {gts}")

        data.groups_data = filtered_genotypes
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
        
        filtered_cavs['value'] = filtered_cavs['logit_es_combined']
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
        
        filtered_cavs['value'] = filtered_cavs['logit_es']
        filtered_cavs['color'] = np.where(filtered_cavs['is_significant'], color, notsignif_color)
        data.cavs_intervals = df_to_variant_intervals(
            filtered_cavs,
            extra_columns=['value', 'color']
        )
        return data


class ReadsParser:
    @staticmethod
    def _convert_cutcounts(reads_dict, allele, interval):
        cuts = np.zeros(len(interval), dtype=np.float32)
        for sample_reads in reads_dict.values():
            cuts += sample_reads[allele]['+']
            cuts += sample_reads[allele]['-']
        return cuts

    @staticmethod
    def _convert_read_list(reads_dict, allele):
        out = []
        for sid, sample_reads in reads_dict.items():
            for read in sample_reads[allele]["fragments"]:
                read.sample_id = sid
                read.base = allele
                out.append(read)
        return out


class ReadsLoader(PlotDataLoader, ReadsParser):
    def _load(self, data: DataBundle, samples_metadata: pd.DataFrame):

        groups_data: pd.Series = data.groups_data
        
        reads = {}
        for sid in groups_data.index:
            cram_path = samples_metadata.loc[sid, "cram_file"]
            reads[sid] = self.get_reads_from_cram(
                cram_path,
                data.interval,
            )
        
        cuts = self._convert_cutcounts(reads, 'N', data.interval)
        reads = self._convert_read_list(reads, 'N')
        data.cutcount_tracks = [{"allele": "N", "cutcounts": cuts}]
        data.reads = np.array(reads)
        return data
    
    def get_reads_from_cram(
        self,
        cram_path,
        interval: GenomicInterval,
    ):
        from footprint_tools.cutcounts import BamFileExtractor
        extractor = BamFileExtractor(cram_path)
        result = {'N': extractor.lookup(interval)}
        extractor.close()
        # {'N': {'+': np.array, '-': np.array, 'fragments': List[GenomicInterval]}}

        return result


class AllelicReadsLoader(PlotDataLoader, ReadsParser):
    def _load(self, data: DataBundle, samples_metadata: pd.DataFrame):

        groups_data: pd.DataFrame = data.groups_data
        interval: GenomicInterval = data.interval

        assert 'parsed_genotype' in groups_data.columns, "groups_data must contain 'parsed_genotype' column"
        gt = np.unique(groups_data['parsed_genotype'])
        assert gt.size == 1, "groups_data must contain only one genotype group"
        variant_interval: VariantInterval = data.variant_interval
        assert variant_interval.overlaps(interval), f"variant_interval must overlap interval. Got {variant_interval.to_str()} and {interval.to_ucsc()}"
        gt = gt[0]

        reads = {}
        for sid in groups_data.index:
            cram_path = samples_metadata.loc[sid, "cram_file"]
            reads[sid] = self.allelic_reads_from_cram(
                cram_path,
                interval,
                variant_interval,
            )

        allele_ref = variant_interval.ref
        allele_alt = variant_interval.alt

        ref_cuts = self._convert_cutcounts(reads, allele_ref, data.interval)
        alt_cuts = self._convert_cutcounts(reads, allele_alt, data.interval)

        if gt == "AA":
            cutcount_tracks = [
                {"allele": allele_ref, "cutcounts": ref_cuts}
            ]
        elif gt == "BB":
            cutcount_tracks = [
                {"allele": allele_alt, "cutcounts": alt_cuts}
            ]
        elif gt == "AB":
            cutcount_tracks = [
                {"allele": allele_ref, "cutcounts": ref_cuts},
                {"allele": allele_alt, "cutcounts": alt_cuts}
            ]
        else:
            raise ValueError(f"Unknown genotype: {gt}")

        data.reads = np.array(
            self._convert_read_list(reads, allele_ref)
            + self._convert_read_list(reads, allele_alt)
        )

        data.cutcount_tracks = cutcount_tracks

        return data
    
    def allelic_reads_from_cram(
        self,
        cram_path,
        interval: GenomicInterval,
        variant_interval: VariantInterval,
    ):
        # Workaround to get the same output as footprint_tools cutcounts extractor
        with AllelicReadsExtractor(cram_path) as extractor:
            reads = extractor[variant_interval] # array of genomic intervals

        result = {
            variant_interval.ref: {
                '+': np.zeros(len(interval), dtype=np.float32),
                '-': np.zeros(len(interval), dtype=np.float32),
                'fragments': []
            },
            variant_interval.alt: {
                '+': np.zeros(len(interval), dtype=np.float32),
                '-': np.zeros(len(interval), dtype=np.float32),
                'fragments': []
            },
        }
        for read_interval in reads:
            allele = read_interval.base
            result[allele]['fragments'].append(read_interval)
            offset_start = read_interval.start - interval.start
            offset_end = read_interval.end - interval.start

            strand = '+' if not read_interval.is_reverse else '-'
            cutcounts_array = result[allele][strand]
            if offset_start >= 0:
                cutcounts_array[offset_start] += 1
            if offset_end < len(interval):
                cutcounts_array[offset_end] += 1
        return result




# DEFUNC
# from genome_tools.data.get_variant_resolved_reads import extract_allelic_reads
# class AllelicReadsLoader(PlotDataLoader):

#     def _load(self, data: DataBundle, samples_metadata: pd.DataFrame, variant_interval: VariantInterval, sample_ids=None):
#         assert variant_interval.overlaps(data.interval), f"variant_interval must overlap data.interval. Got {variant_interval.to_str()} and {data.interval.to_ucsc()}"

#         variant_data: pd.DataFrame = data.variant_genotypes.query('parsed_genotype == "AB"')
#         variant_data = filter_df_to_interval(variant_data, variant_interval)
#         if variant_data.empty:
#             raise ValueError("No heterozygous genotypes found for the specified variant_interval.")
#         assert variant_data.index.is_unique

#         if sample_ids is None:
#             sample_ids = variant_data.index.tolist()
#         elif isinstance(sample_ids, (str, int, float)):
#             sample_ids = [sample_ids]

#         cram_paths = samples_metadata.loc[sample_ids, 'cram_file']

#         reads = {}
#         for sample_id, cram_path in tqdm(zip(sample_ids, cram_paths), total=len(sample_ids), desc="Allelic reads loader"):
#             extracted_reads = extract_allelic_reads(cram_path, variant_interval, data.interval)
#             if not self.check_reads(extracted_reads, variant_interval):
#                 continue
#             reads[sample_id] = extracted_reads
#         data.reads = reads
#         return data
    
#     def check_reads(self, extracted_reads: List[GenomicInterval], variant_interval: VariantInterval):
#         n_ref, n_alt = 0, 0
#         for read in extracted_reads:
#             n_ref += read.base == variant_interval.ref
#             n_alt += read.base == variant_interval.alt
#         if n_ref > 0 and n_alt > 0:
#             return True
#         return False