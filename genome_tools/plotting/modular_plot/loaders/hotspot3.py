import pandas as pd


from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle

from genome_tools.data.extractors import ChromParquetExtractor

from hotspot3.io.readers import ChromReader
from hotspot3.background_fit.fit import GlobalBackgroundFit


class AggCutcountsLoader(PlotDataLoader):
    def _load(self, data: DataBundle, cutcounts_tabix_file, mappable_bases_tabix, chrom_sizes: pd.DataFrame):
        reader = ChromReader(chrom_sizes.loc[data.interval.chrom])
        signal = reader.extract_mappable_agg_cutcounts(
            cutcounts_file=cutcounts_tabix_file,
            mappable_file=mappable_bases_tabix
        )
        data.signal = signal
        return data


class PerBpBackgroundTrackLoader(PlotDataLoader):
    def _load(self, data: DataBundle, parquet_bg_file):
        with ChromParquetExtractor(parquet_bg_file, columns=['sliding_r', 'sliding_p', 'tr']) as pqt:
            df: pd.DataFrame = pqt[data.interval]
            data.mean_bg = df.eval('sliding_r * sliding_p / (1 - sliding_p)').values
            data.fit_threshold = df['tr'].values
        return data


class HighSignalMaskLoader(PlotDataLoader):
    def _load(self, data: DataBundle):
        g_fit = GlobalBackgroundFit()
        data.high_signal_mask = g_fit.get_signal_mask_for_tr(
            data.signal, 
            data.fit_threshold,
            75
        )
        return data
