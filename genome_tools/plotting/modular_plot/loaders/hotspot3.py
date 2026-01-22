import pandas as pd
import numpy as np

from scipy.signal import find_peaks

from hotspot3.io.readers import ChromReader
from hotspot3.background_fit.fit import GlobalBackgroundFit

from genome_tools import GenomicInterval
from genome_tools.data.extractors import ChromParquetExtractor, TabixExtractor


from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle

from genome_tools.plotting.modular_plot.loaders.basic import SignalLoader


class AggCutcountsLoader(PlotDataLoader):
    def _load(self, data: DataBundle, cutcounts_tabix_file, mappable_bases_tabix, chrom_sizes: pd.DataFrame):
        chrom_interval = GenomicInterval(
            data.interval.chrom,
            0,
            chrom_sizes.loc[data.interval.chrom, 'len']
        )
        reader = ChromReader(chrom_interval)
        signal = reader.extract_mappable_agg_cutcounts(
            cutcounts_file=cutcounts_tabix_file,
            mappable_file=mappable_bases_tabix
        )[data.interval.start:data.interval.end]
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


class SmoothedSignalLoader(PlotDataLoader):
    def _load(self, data: DataBundle, smoothed_signal_parquet):
        with ChromParquetExtractor(smoothed_signal_parquet, columns=['smoothed']) as pqt:
            data.smoothed_signal =  pqt[data.interval]['smoothed'].values
        data.maxima = find_peaks(data.smoothed_signal)[0]
        return data


class BackgroundDensityLoader(PlotDataLoader):
    def _load(self, data: DataBundle, bg_tabix, read_depth):
        with TabixExtractor(bg_tabix) as extractor:
            bg_df = extractor[data.interval]
            bg_df['mean'] = bg_df.eval('bg_r * bg_p / (1 - bg_p)').values # counts
            bg_df['mean_density'] = bg_df['mean'] * 1_000_000 / read_depth # density
        L = data.interval.end - data.interval.start
        bg = np.zeros(L, dtype=np.float32)
        starts = (bg_df['start'].to_numpy() - data.interval.start).clip(0, L)
        ends   = (bg_df['end'].to_numpy() - data.interval.start).clip(0, L)
        vals   = bg_df['mean_density'].to_numpy()

        for s, e, v in zip(starts, ends, vals):
            bg[s:e] = v

        data.background_density = bg
        return data
