from bx.align import maf

import numpy as np
from tqdm import tqdm
from genome_tools import GenomicInterval
from genome_tools.data.utils import remap_columns
import pandas as pd


class BetweenSpeciesMap:
    def __init__(self, mapping, root_interval: GenomicInterval, root_species='Species1', target_species='Species2'):
        """
        General init of the class. Usually, not called directly. Use .from_maf class method to instantiate the class instead.
        mapping - mapping dict in the following format {root_chrom: {root_pos: (target_chrom, target_pos)}}. Positions are 0-based. Usually generated from .maf file (see .from_maf method)
        root_interval - corresponding genomic region GenomicInterval
        """
        self.mapping = mapping

        self.reverse = self._build_reverse(mapping)

        self.root_interval = root_interval

        self.root_species = root_species
        self.target_species = target_species

    def __repr__(self):
        return f'BetweenSpeciesMap({self.root_species} -> {self.target_species}, interval={self.root_interval.to_ucsc()})'

    @staticmethod
    def _build_reverse(forward):
        rev = {}
        for r_chr, pos_dict in forward.items():
            for r_pos, (t_chr, t_pos) in pos_dict.items():
                rev.setdefault(t_chr, {})[t_pos] = (r_chr, r_pos)
        return rev

    @classmethod
    def from_maf(cls, maf_path, root_interval, root_species='Homo_sapiens', target_species='Mus_musculus'):
        """
        maf_path - result of hal2maf for a region
        """
        mapping = {}

        if target_species == root_species:
            positions = np.arange(root_interval.start, root_interval.end)
            positions = {x: (root_interval.chrom, x) for x in positions}
            mapping = {root_interval.chrom: positions}

        else:
            skipped_all = True
            with open(maf_path) as f:
                reader = maf.Reader(f)

                for block in reader:
                    root = None
                    targets = []

                    for comp in block.components:
                        if comp.src.startswith(root_species + '.'):
                            root = comp
                        elif comp.src.startswith(target_species + '.'):
                            targets.append(comp)
                            skipped_all = False

                    if not root or not targets:
                        continue

                    # Prefer the longest component; use shorter components
                    # only for root positions not already mapped.
                    targets = sorted(targets, key=lambda x: x.size, reverse=True)

                    for target in targets:

                        r_chrom = root.src.removeprefix(root_species + '.')
                        t_chrom = target.src.removeprefix(target_species + '.')

                        r_seq = root.text
                        t_seq = target.text

                        if root.strand == '+':
                            r_pos = root.start
                            r_step = 1
                        else:
                            r_pos = root.src_size - root.start - 1
                            r_step = -1

                        if target.strand == '+':
                            t_pos = target.start
                            t_step = 1
                        else:
                            t_pos = target.src_size - target.start - 1
                            t_step = -1

                        for i in range(len(r_seq)):
                            cur_r = None
                            cur_t = None

                            if r_seq[i] != '-':
                                cur_r = r_pos
                                r_pos += r_step

                            if t_seq[i] != '-':
                                cur_t = t_pos
                                t_pos += t_step

                            if cur_r is not None and cur_t is not None:
                                assert GenomicInterval(r_chrom, cur_r, cur_r + 1).overlaps(root_interval), f'MAF file mapping contains position outside of root_interval {root_interval.to_ucsc()}. Are you sure the interval corresponds to provided MAF file?'
                                mapping.setdefault(r_chrom, {}).setdefault(cur_r, (t_chrom, cur_t))

            if skipped_all:
                raise ValueError(f'Species {target_species} not present in the mapping. Check the spelling.')
        return cls(mapping, root_interval, root_species=root_species, target_species=target_species)

    def map_position_root_to_target(self, chrom, pos):
        return self.mapping.get(chrom, {}).get(pos)

    def map_position_target_to_root(self, chrom, pos):
        return self.reverse.get(chrom, {}).get(pos)


    def map_interval_to_root(self, interval: GenomicInterval):
        return self._map_interval(
            interval=interval,
            mapper_method=self.map_position_target_to_root
        )


    def map_interval_to_target(self, interval: GenomicInterval):
        return self._map_interval(
            interval=interval,
            mapper_method=self.map_position_root_to_target
        )

    def map_row_target_to_root(self, row):
        start_res = self.map_position_target_to_root(row['#chr'], row['start'])
        if start_res is not None:
            new_chrom, new_start = start_res
        else:
            new_start = pd.NA
            new_chrom = pd.NA

        end_res = self.map_position_target_to_root(row['#chr'], row['end'] - 1)
        if end_res is not None:
            _, new_end = end_res
            new_end += 1
        else:
            new_end = pd.NA
        if pd.isna(new_start) or pd.isna(new_end) or (new_end - new_start != row['end'] - row['start']):
            new_chrom = pd.NA
            new_start = pd.NA
            new_end = pd.NA

        return pd.Series(
            {
                '#chr': new_chrom,
                'start': new_start,
                'end': new_end
            }
        )

    def map_target_df_to_root(self, df):
        return df.progress_apply(
            self.map_row_target_to_root, axis=1
        )

    def _map_interval(self, interval: GenomicInterval, mapper_method):
        """
        interval: interval to map
        mapper_method: self.map_position_target_to_root or self.map_position_root_to_target
        Returns:
            GenomicInterval or None
        """
        start = interval.start
        end = interval.end
        chrom = interval.chrom

        while start < end:
            m_start = mapper_method(chrom, start)
            m_end = mapper_method(chrom, end - 1)  # correct for half-open

            if m_start is not None and m_end is not None:
                m_chrom_s, m_start = m_start
                m_chrom_e, m_end = m_end

                if m_chrom_s == m_chrom_e:
                    return GenomicInterval(m_chrom_s, min(m_start, m_end), max(m_start, m_end) + 1,
                                           strand='+' if m_start <= m_end else '-')
                else:
                    raise ValueError("Chrom mismatch")

            # shrink only failing side(s)
            if m_start is None:
                print(f"Warning: shrinking start {chrom}:{start}-{end}")
                start += 1

            if m_end is None:
                print(f"Warning: shrinking end {chrom}:{start}-{end}")
                end -= 1

        return None


class ParsedFastaHandler:

    def __init__(
            self,
            parsed_fasta,
            root_target_mapping: BetweenSpeciesMap,
        ):
        self.parsed_fasta = parsed_fasta
        self.root_target_mapping = root_target_mapping
        self.root_interval = root_target_mapping.root_interval

        root_seq = parsed_fasta[root_target_mapping.root_species]

        self.root_pos_to_col = {
            pos: col
            for pos, col in zip(
                range(self.root_interval.start, self.root_interval.end),
                (i for i, b in enumerate(root_seq) if b != '-')
            )
        }
        self.root_pos_to_col[self.root_interval.end] = len(root_seq)

    @classmethod
    def from_fasta(cls, fasta_path, root_target_mapping: BetweenSpeciesMap):
        parsed_fasta = {}
        key = None

        with open(fasta_path) as f:
            for line in tqdm(f):
                line = line.strip()
                if line.startswith('>'):
                    if key is not None:
                        parsed_fasta[key] = sequence
                    key = line[1:].strip()
                    sequence = ""
                else:
                    sequence += line

        if key is not None:
            parsed_fasta[key] = sequence

        return cls(parsed_fasta, root_target_mapping)

    def __getitem__(self, interval):
        if (
            interval.chrom != self.root_interval.chrom
            or interval.start < self.root_interval.start
            or interval.end > self.root_interval.end
        ):
            raise ValueError(
                f'{interval.to_ucsc()} does not fully overlap with '
                f'{self.root_interval.to_ucsc()}'
            )

        start = self.root_pos_to_col[interval.start]
        end = self.root_pos_to_col[interval.end]

        return {
            species: seq[start:end]
            for species, seq in self.parsed_fasta.items()
        }

def map_matrix_to_interval(
        matrix: np.ndarray,
        matrix_interval: GenomicInterval,
        target_interval: GenomicInterval,
        mapping: BetweenSpeciesMap,
    ) -> np.ndarray:
    """Cross-species: source column comes from the position mapping."""
    assert matrix.shape[-1] == len(matrix_interval)
    src_cols = np.full(len(target_interval), -1, dtype=np.intp)
    for i, pos in enumerate(range(target_interval.start, target_interval.end)):
        mapped = mapping.map_position_root_to_target(target_interval.chrom, pos)
        if mapped is not None:
            src_cols[i] = mapped[1] - matrix_interval.start

    if not (src_cols >= 0).any():
        print(f"Warning: mapping produced no valid columns for {target_interval} with {matrix_interval}. Check the mapping and interval overlap.")

    return remap_columns(matrix, src_cols)
