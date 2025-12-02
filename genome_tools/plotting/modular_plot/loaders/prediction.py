import numpy as np
import anndata as ad

import torch
from torch.utils.data import DataLoader

from vinson.datasets.sequence import SequenceEmbedDataset
from vinson.utils.data_formatting import extract_data_from_backed_anndata
from vinson.postprocessing.interpretation import ModelWrapper

from genome_tools import GenomicInterval

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle


class AttributionsLoader(PlotDataLoader):
    def _load(self, data: DataBundle,
                sample_id: str,
                dhs_id: str,
                anndata: ad.AnnData,
                model_config: dict,
                model_wrapper: ModelWrapper,
                fasta_file: str,
                genotype_file: str,
                print_convergence_deltas=False,
        ):
        X, X_embed, interval = self.make_data_for_model(
            sample_id=sample_id,
            dhs_id=dhs_id,
            anndata=anndata,
            model_config=model_config,
            fasta_file=fasta_file,
            genotype_file=genotype_file,
        )
        assert interval.overlaps(data.interval)

        attrs = model_wrapper.get_sequence_attributions(X, X_embed, print_convergence_deltas=print_convergence_deltas).squeeze(0).numpy()

        data.matrix = self.align_matrix_to_interval(
            matrix=attrs,
            matrix_interval=interval,
            target_interval=data.interval,
        ).T
        data.sequence_weights = data.matrix.sum(axis=1)
        return data
    
    @staticmethod
    def make_data_for_model(
        sample_id,
        dhs_id,
        anndata,
        model_config: dict,
        fasta_file: str,
        genotype_file: str,
    ):
        dataset_kwargs: dict = model_config['data_params']
        dataset_kwargs.update(
            dict(
                reverse_complement=False,
                jitter=0,
                noise=0,
            )
        )
        
        data = extract_data_from_backed_anndata(
            anndata,
            sample_ids=[sample_id],
            dhs_ids=[dhs_id],
        )


        dataset = SequenceEmbedDataset(
            data=data,
            embeddings_df=anndata[[sample_id], [dhs_id]].obsm['motif_embeddings'],
            fasta_file=fasta_file,
            genotype_file=genotype_file,
            **dataset_kwargs,
        )
        
        dataloader = DataLoader(
            dataset,
            batch_size=len(dataset),
            shuffle=False,
            num_workers=8,
            pin_memory=True if torch.cuda.is_available() else False,
            drop_last=False,
        )

        batch = next(iter(dataloader))

        summit, chrom = anndata.var.loc[dhs_id, ['dhs_summit', '#chr']]
        interval = GenomicInterval(chrom, summit, summit, dhs_summit=summit).widen(672)

        return batch['ohe_seq'], batch['embed'], interval
    
    @staticmethod
    def align_matrix_to_interval(
        matrix: np.ndarray,
        matrix_interval: GenomicInterval,
        target_interval: GenomicInterval,
    ):
        assert matrix.shape[1] == len(matrix_interval)
        aligned = np.zeros(
            (matrix.shape[0], len(target_interval)),
            dtype=matrix.dtype
        )
        offset = matrix_interval.start - target_interval.start
        t_start = max(offset, 0)
        t_end = min(offset + len(matrix_interval), len(target_interval))

        if t_start >= t_end:
            return aligned  # no overlap

        m_start = max(-offset, 0)
        m_end = m_start + (t_end - t_start)

        aligned[:, t_start:t_end] = matrix[:, m_start:m_end]
        return aligned
