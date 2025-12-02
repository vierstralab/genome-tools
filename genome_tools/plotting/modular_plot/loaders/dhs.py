import pandas as pd

from genome_tools.utils.signal import smooth_and_aggregate_per_nucleotide_signal

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle
from .basic import SegmentsLoader


class DHSIndexLoader(SegmentsLoader):
    __intervals_attr__ = 'dhs_intervals'
    
    def _load(self, data: DataBundle, dhs_index: pd.DataFrame, extra_columns=None, rectprops_columns=None):
        return super()._load(
            data, 
            segments_df=dhs_index, 
            extra_columns=extra_columns,
            rectprops_columns=rectprops_columns
        )


class ComponentTracksLoader(PlotDataLoader):

    def _load(self, data: DataBundle, cutcounts_files, smooth=True, step=20, bandwidth=150, nmf_components=None):
        if not smooth:
            bandwidth = 1

        if nmf_components is None:
            nmf_components = cutcounts_files.keys()

        data.nmf_components = nmf_components

        data.component_tracks = []
        for component in nmf_components:
            segs = smooth_and_aggregate_per_nucleotide_signal(data.interval,
                                                              cutcounts_files[component],
                                                              step=step, bandwidth=bandwidth)
            data.component_tracks.append(segs)
        return data


class DHSLoadingsLoader(PlotDataLoader):
    required_loader_kwargs = ['H']
