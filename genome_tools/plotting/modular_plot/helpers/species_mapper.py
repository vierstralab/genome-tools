from bx.align import maf
import os
from tqdm import tqdm
from genome_tools import GenomicInterval
import pandas as pd


class BetweenSpeciesMap:

    def __init__(self, mapping):
        """
        E.g. {human_chrom: {human_pos: (mouse_chrom, mouse_pos)}}
        """
        self.mapping = mapping

        self.reverse = self._build_reverse(mapping)

    @staticmethod
    def _build_reverse(forward):
        rev = {}
        for r_chr, pos_dict in forward.items():
            for r_pos, (t_chr, t_pos) in pos_dict.items():
                rev.setdefault(t_chr, {})[t_pos] = (r_chr, r_pos)
        return rev

    @classmethod
    def from_maf(cls, maf_path, root_species='Homo_sapiens', target_species='Mus_musculus'):
        mapping = {}

        with open(maf_path) as f:
            reader = maf.Reader(f)

            for block in reader:
                human = None
                mouse = None

                for comp in block.components:
                    if comp.src.startswith(root_species):
                        human = comp
                    elif comp.src.startswith(target_species):
                        mouse = comp

                if not (human and mouse):
                    continue

                h_chrom = human.src.split('.')[-1]
                m_chrom = mouse.src.split('.')[-1]

                h_seq = human.text
                m_seq = mouse.text

                h_pos = human.start
                m_pos = mouse.start

                h_step = 1 if human.strand == '+' else -1
                m_step = 1 if mouse.strand == '+' else -1

                for i in range(len(h_seq)):
                    cur_h = None
                    cur_m = None

                    if h_seq[i] != '-':
                        cur_h = h_pos
                        h_pos += h_step

                    if m_seq[i] != '-':
                        cur_m = m_pos
                        m_pos += m_step

                    if cur_h is not None and cur_m is not None:
                        mapping.setdefault(h_chrom, {})[cur_h] = (m_chrom, cur_m)

        return cls(mapping)

    def map_pos_root_to_target(self, chrom, pos):
        return self.mapping.get(chrom, {}).get(pos)

    def map_pos_target_to_root(self, chrom, pos):
        return self.reverse.get(chrom, {}).get(pos)
    
    def map_interval_to_root(self, interval: GenomicInterval):
        return self.map_pos_interval(
            interval=interval,
            mapper_method=self.map_pos_target_to_root
        )
    
    def map_interval_to_target(self, interval: GenomicInterval):
        return self.map_pos_interval(
            interval=interval,
            mapper_method=self.map_pos_root_to_target
        )

    def map_row(self, row):
        start_res = self.map_pos_target_to_root(row['#chr'], row['start'])
        if start_res is not None:
            new_chrom, new_start = start_res
        else:
            new_start = pd.NA
            new_chrom = pd.NA
        
        end_res = self.map_pos_target_to_root(row['#chr'], row['end'] - 1)
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
            self.map_row, axis=1
        )

    def map_pos_interval(self, interval: GenomicInterval, mapper_method):
        """
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
                    return GenomicInterval(m_chrom_s, m_start, m_end + 1)
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
            human_mouse_mapping: BetweenSpeciesMap,
            initial_region: GenomicInterval
        ):
        self.parsed_fasta = parsed_fasta
        self.initial_region = initial_region
        self.human_mouse_mapping = human_mouse_mapping

    @classmethod
    def from_folder(cls, fpath, initial_region, human_mouse_mapping=None):
        if os.path.isdir(fpath):
            fname = f"{fpath}/{initial_region.to_ucsc()}.fasta"
            human_mouse_mapping = f"{fpath}/{initial_region.to_ucsc()}.maf"
        else:
            fname = fpath
            assert human_mouse_mapping is not None
        parsed_fasta = {}
        key = None
        
        
        with open(fname) as f:
            for line in tqdm(f.readlines()):
                line = line.strip()
                if line.startswith('>'):
                    if key is not None:
                        parsed_fasta[key] = sequence
                    key = line[2:]
                    sequence = ""
                else:
                    sequence = sequence + line
        
        if not isinstance(human_mouse_mapping, BetweenSpeciesMap):
            human_mouse_mapping = BetweenSpeciesMap.from_maf(human_mouse_mapping)
        return cls(parsed_fasta, human_mouse_mapping, initial_region)

    def __getitem__(self, interval):
        has_overlap = (self.initial_region.start <= interval.start) and (self.initial_region.end >= interval.end)
        if not interval.overlaps(self.initial_region) or not has_overlap:
            raise ValueError(f'{interval.to_ucsc()} does not fully overlap with {self.initial_region.to_ucsc()}')
        offset = interval.start - self.initial_region.start
        length = interval.end - interval.start

        return {
            x: y[offset:offset + length]
            for x, y
            in self.parsed_fasta.items()
        }
