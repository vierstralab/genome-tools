import numpy as np
import pandas as pd

from collections import defaultdict
from genome_tools.data.extractors import FastaExtractor
from genome_tools import GenomicInterval, df_to_genomic_intervals, genomic_intervals_to_df


from typing import Union, List

import pickle
from tqdm import tqdm

try:
    from hotspot3.helpers.utils import correct_offset
    import bottleneck as bn
except ImportError:
    print('Please install hotspot3 to use sampling')


def sanitize_bed_data(data: Union[List[GenomicInterval], pd.DataFrame]):
    if isinstance(data, pd.DataFrame):
        data = df_to_genomic_intervals(data)
    return data


class GCTrack:
    def __init__(self, chrom_sizes: pd.DataFrame, fasta_file, window_size=1345, n_bins=100):
        self.chrom_sizes = chrom_sizes
        self.chrom_names = chrom_sizes.index.tolist()
        self.fasta_file = fasta_file

        self.window_size = window_size
        self.flank_length = (window_size - 1) // 2

        self.n_bins = n_bins
        self.gc_bins = np.linspace(0, 1, n_bins + 1)

        self.track = {
            chrom: np.full((length, 2), -1, dtype=np.int16)
            for chrom, length in zip(self.chrom_names, self.chrom_sizes['length'])
        }
        self.sampling_index = defaultdict(pd.DataFrame)

        self.built = False
        self.masked = False
        self.has_sampling_index = False
    

    def parse_key(self, key):
        assert isinstance(key, tuple)
        if len(key) != 2:
            raise TypeError("Expected (GenomicInterval, mode)")
        interval, mode = key
  

        if not isinstance(interval, GenomicInterval):
            raise TypeError("First argument must be a GenomicInterval")
        if mode not in ('masked', 'raw', 'start', 'end'):
            raise ValueError("Mode must be 'masked' or 'raw'")
        
        if mode == 'masked':
            col = 1
        elif mode == 'raw':
            col = 0
        else:
            raise ValueError("Invalid mode")
        return interval, col

    def __getitem__(self, key):
        interval, col = self.parse_key(key)
        return self.track[interval.chrom][interval.start:interval.end, col]
    
    def __setitem__(self, key, value):
        interval, col = self.parse_key(key)
        self.track[interval.chrom][interval.start:interval.end, col] = value
    
    def at(self, key):
        assert isinstance(key, GenomicInterval), "Key must be a GenomicInterval"
        center = key.center
        return self[center, 'raw'][0]
    
    @correct_offset
    def bn_mean_correct_offset(self, array, window, *args, **kwargs):
        return bn.move_mean(array, window, *args, **kwargs).astype(np.float32)
    

    def build(self):
        extractor = FastaExtractor(self.fasta_file)
        for chrom, row in self.chrom_sizes.iterrows():
            length = int(row['length'])
            gi = GenomicInterval(chrom, 0, length)
            seq = extractor[gi].upper()
            arr = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)

            is_gc = np.select(
                [
                    (arr == ord('G')) | (arr == ord('C')),
                    (arr == ord('A')) | (arr == ord('T')),
                ],
                [
                    1,
                    0
                ],
                default=np.nan  # Other bases
            ).astype(np.float32)


            # GC fraction per window
            print('Convolving')
            print(self.window_size)
            gc_frac = self.bn_mean_correct_offset(is_gc, self.window_size, self.flank_length).astype(np.float32)
            bins = np.digitize(gc_frac, self.gc_bins[1:-1])
            print('Convolution done')
            self[gi, 'raw'] = bins
            self[gi, 'masked'] = bins
        extractor.close()

        self.built = True

    def mask(self, exclude_regions: Union[List[GenomicInterval], pd.DataFrame]):
        # Group excluded intervals by chromosome
        grouped = defaultdict(list)
        for ex in sanitize_bed_data(exclude_regions):
            grouped[ex.chrom].append(ex)

        for chrom, chrom_track in self.track.items():
            mask = np.zeros(len(chrom_track), dtype=np.float32)
            ## TODO: double check
            for ex in grouped.get(chrom, []):
                # Mark the base-level exclusion mask
                 # ensure in-bounds 
                start = max(0, ex.start)
                end = min(ex.end, start + len(mask)) 
                mask[start:end] = 1

            # Convolve once per chromosome
            exclude_windows = (self.bn_mean_correct_offset(
                mask,
                self.window_size
            ) > 0).astype(bool)
            chrom_track[:, 1] = chrom_track[:, 0]
            chrom_track[exclude_windows > 0, 1] = -1
        self.masked = True


    def build_sampling_index(self):
        if not self.masked:
            print('Warning! No regions were masked. All data will be used in sampling.')
        sampling_data = defaultdict(list)
        for chromosome in self.chrom_names:
            gc_bins_track = self.track[chromosome][:, 1]
            sorted_idx = np.argsort(gc_bins_track)
            sorted_vals = gc_bins_track[sorted_idx]
            unique_gc_bins, start_idx = np.unique(sorted_vals, return_index=True)

            for gc_bin, start, start_next in zip(
                unique_gc_bins,
                start_idx,
                np.append(start_idx[1:], len(gc_bins_track))
            ):
                if pd.isna(gc_bin) or gc_bin == -1:
                    continue
                data = pd.DataFrame(
                    dict(center=sorted_idx[start:start_next], chrom=chromosome)
                )
                sampling_data[gc_bin].append(data)
        
        for key, value in sampling_data.items():
            self.sampling_index[key] = pd.concat(value, ignore_index=True)
        self.has_sampling_index = True


class GCSampler:
    def __init__(self, gc_track: GCTrack, seed: int = 777, max_iter: int = 1000):
        self.gc_track = gc_track
        self.chrom_sizes = gc_track.chrom_sizes
        self.window_size = gc_track.window_size
        self.flank_length = self.gc_track.flank_length

        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.max_iter = max_iter

    def sample(self, positives: Union[List[GenomicInterval], pd.DataFrame], n=1, *args, **kwargs) -> pd.DataFrame:
        if self.gc_track.built is False:
            raise ValueError("GC track is not built. Call `build` method for GC track first.")
        
        if not self.gc_track.masked:
            print('Warning! No regions are masked in GC track. All data will be used in sampling.')
        
        if not self.gc_track.has_sampling_index:
            print('GC track doesn\'t have a sampling index. Call `build_sampling_index` for GC track to enable efficent sampling')
            return self.rejection_sampling(positives, n, *args, **kwargs)
        else:
            return self.from_index_sampling(positives, n, *args, **kwargs)

    def from_index_sampling(self, positives: list[GenomicInterval], n: int):
        gc_bins, gc_bin_counts = np.unique(
            [self.gc_track.at(positive) for positive in positives],
            return_counts=True
        )
        gc_bin_counts = gc_bin_counts * n

        result = []
        for bin, count in zip(gc_bins, gc_bin_counts):
            if len(self.gc_track.sampling_index[bin]) < count:
                print(f"Warning: Not enough samples for GC bin {bin}. Requested {count}, available {len(self.gc_track.sampling_index[bin])}. Sampling with replacement")
                replace = True
            else:
                replace = False
            result.append(
                self.gc_track.sampling_index[bin].sample(
                    n=count,
                    replace=replace,
                    random_state=self.seed
                )
            )
        result = pd.concat(result, ignore_index=True)
        result['start'] = result['center'] - self.flank_length

        result['end'] = result['center'] + self.flank_length + 1
        return result

    def rejection_sampling(
            self,
            positives: list[GenomicInterval],
            n: int, 
            max_iter: int=1000,
            same_chromosome=False
        ) -> pd.DataFrame:
        """
        Sample multiple genomic intervals matching the GC track.
        """
        result = []
        for positive in tqdm(positives):
            gc_bin = self.gc_track.at(positive)
            for _ in range(n):
                interval = self._rejection_sample_one(
                    positive,
                    gc_bin,
                    max_iter=max_iter,
                    same_chromosome=same_chromosome
                )
                if interval not in result and interval is not None:
                    result.append(interval)
            
        return genomic_intervals_to_df(result).drop(columns='name')

    def _rejection_sample_one(
            self,
            reference_interval: GenomicInterval,
            bin_id: int,
            max_iter: int,
            same_chromosome: bool
        ) -> int:
        """
        Sample a single genomic interval matching the GC track.
        """
        for _ in range(max_iter):
            if same_chromosome:
                chrom = reference_interval.chrom
            else:
                chrom = self.rng.choice(self.chrom_sizes.index)
            sampled_center = self.rng.integers(
                self.flank_length,
                self.chrom_sizes.at[chrom, 'length'] - self.flank_length
            )
            interval = GenomicInterval(chrom, sampled_center, sampled_center + 1)
            if self.gc_track[interval, 'masked'][0] == bin_id:
                return interval.widen(self.flank_length)
        else:
            print(f"Warning: Maximum iterations {max_iter} reached without finding a valid sample.")
