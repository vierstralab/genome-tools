from typing import List
import pandas as pd
import numpy as np

from genome_tools.data.extractors import TabixExtractor, BigwigExtractor
from genome_tools.utils.signal import smooth_and_aggregate_per_nucleotide_signal
from genome_tools import GenomicInterval, df_to_genomic_intervals

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class AnnotationRegionsLoader(PlotDataLoader):
    def _load(self, data: DataBundle, annotation_regions: List[GenomicInterval]):
        for region in annotation_regions:
            assert region.overlaps(data.interval), f"Annotation regions must overlap with interval. Got non-overlapping region {region.to_ucsc()} (interval={data.interval.to_ucsc()})"

        data.annotation_regions = annotation_regions
        return data


# TODO: wrap as class factory
class IdeogramLoader(PlotDataLoader):
    required_loader_kwargs = ['ideogram_data']


class GencodeLoader(PlotDataLoader):
    required_loader_kwargs = ['gencode_annotation_file']


class SignalLoader(PlotDataLoader):
    def _load(self, data, signal_file):
        with BigwigExtractor(signal_file) as extractor:
            signal = np.nan_to_num(extractor[data.interval])
        data.signal = signal
        return data

class AverageSignalLoader(PlotDataLoader):
    
    def _load(self, data, signal_files, smooth=True, step=20, bandwidth=150):
        if not smooth:
            bandwidth = 1
        segs = smooth_and_aggregate_per_nucleotide_signal(
            data.interval,
            signal_files,
            step=step,
            bandwidth=bandwidth
        )
        data.signal = segs
        return data


class SegmentsLoader(PlotDataLoader):
    __intervals_attr__ = 'intervals'

    def _load(self, data: DataBundle, segments_df: pd.DataFrame, extra_columns=None, rectprops_columns=None):
        if rectprops_columns is None:
            rectprops_columns = []
        if extra_columns is None:
            extra_columns = []
        setattr(
            data,
            self.__intervals_attr__, 
            df_to_genomic_intervals(
                segments_df.reset_index(drop=True).reset_index(),
                data.interval,
                extra_columns=['index'] + extra_columns + rectprops_columns
            )
        )

        if rectprops_columns:
            for interval in getattr(data, self.__intervals_attr__):
                interval.rectprops = {
                    col: getattr(interval, col) for col in rectprops_columns
                    }
        return data


class SegmentsTabixLoader(SegmentsLoader):
    def _load(self, data: DataBundle, tabix_file: str, extra_columns=None, rectprops_columns=None, **tabix_loader_kwargs):
        with TabixExtractor(tabix_file, **tabix_loader_kwargs) as extractor:
            segments_df = extractor[data.interval]
        return super()._load(
            data,
            segments_df=segments_df,
            extra_columns=extra_columns,
            rectprops_columns=rectprops_columns
        )

