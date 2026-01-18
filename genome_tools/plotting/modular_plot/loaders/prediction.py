import numpy as np
import anndata as ad

from scipy.interpolate import interp1d

from typing import Tuple

import torch
from torch.utils.data import DataLoader

from vinson.datasets.sequence import SequenceEmbedDataset
from vinson.utils.data_formatting import extract_data_from_backed_anndata, VinsonData
from vinson.postprocessing.interpretation import ModelWrapper

from genome_tools import GenomicInterval

from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle



class PredictionDataLoader(PlotDataLoader):
    @staticmethod
    def from_backed_anndata(
        anndata: ad.AnnData,
        sample_id: str,
        dhs_id: str,
    ) -> Tuple[VinsonData, GenomicInterval]:
        data = extract_data_from_backed_anndata(
            anndata,
            sample_ids=[sample_id],
            dhs_ids=[dhs_id],
        )

        summit, chrom = anndata.var.loc[dhs_id, ['dhs_summit', '#chr']]
        interval = GenomicInterval(chrom, summit, summit, dhs_summit=summit).widen(672)
        return data, interval

    @staticmethod
    def from_prediction_coordinates(
        anndata: ad.AnnData,
        sample_id: str,
        coordinates: np.ndarray,
        chrom: str,
    ):
        
        raw_data = {
            'sample_id': np.array([sample_id] * len(coordinates), dtype=np.str_),
            'dhs_id': np.array([''] * len(coordinates), dtype=np.str_),
            'read_depth': np.array([anndata.obs.loc[sample_id, 'nuclear_reads']] * len(coordinates)),
            'chrom': np.array([chrom] * len(coordinates), dtype=np.str_),
            'summit': np.array(coordinates),
            'background': np.array([0.0] * len(coordinates), dtype=np.float32), # dummy
            'class': np.array([0] * len(coordinates), dtype=np.int8), # dummy
            'density': np.array([0.0] * len(coordinates)), # dummy
        }

        return VinsonData.from_raw(
            raw_data=raw_data,
            is_variant=False,
            embeddings_df=anndata[[sample_id], :].obsm['motif_embeddings'],
        )

    @staticmethod
    def get_dataset(
        data: VinsonData,
        model_config: dict,
        fasta_file: str,
        genotype_file: str=None,
    ):
        # FIXME: avoid reading the whole dataset
        dataset_kwargs: dict = model_config['data_params']
        dataset_kwargs.update(
            dict(
                reverse_complement=False,
                jitter=0,
                noise=0,
            )
        )

        dataset = SequenceEmbedDataset(
            data=data,
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
        return batch


class PredictedSignalLoader(PredictionDataLoader):
    def _load(self, data: DataBundle,
                sample_id: str,
                anndata: ad.AnnData,
                model_config: dict,
                model_wrapper: ModelWrapper,
                fasta_file: str,
                genotype_file: str=None,
                interp1d_kind='linear',
                step=20,
        ):
        initial_interval: GenomicInterval = data.interval

        prediction_starts = np.arange(
            initial_interval.start,
            initial_interval.end + step,
            step
        )

        model_data = self.from_prediction_coordinates(
            anndata,
            sample_id=sample_id,
            coordinates=prediction_starts,
            chrom=initial_interval.chrom,
        )

        batch = self.get_dataset(
            model_data,
            model_config,
            fasta_file,
            genotype_file,
        )
        X_seq = batch["ohe_seq"].to(model_wrapper.model.device, non_blocking=True)
        X_embed = batch["embed"].to(model_wrapper.model.device, non_blocking=True)

        pred_density = model_wrapper(X_seq, X_embed).detach().cpu().numpy().squeeze()

        full_positions = np.arange(initial_interval.start, initial_interval.end)

        data.signal = interp1d(
            prediction_starts,
            pred_density,
            kind=interp1d_kind,
            assume_sorted=True
        )(full_positions)
        return data


class AttributionsLoader(PredictionDataLoader):
    def _load(self, data: DataBundle,
                sample_id: str,
                dhs_id: str,
                anndata: ad.AnnData,
                model_config: dict,
                model_wrapper: ModelWrapper,
                fasta_file: str,
                genotype_file: str,
                print_convergence_deltas=False,
                n_shuffles=20,
        ):
        input_data, interval = self.from_backed_anndata(
            anndata,
            sample_id=sample_id,
            dhs_id=dhs_id,
        )
        assert interval.overlaps(data.interval), f"Data interval {data.interval} does not overlap with the dhs interval {interval}, {dhs_id}"

        batch = self.get_dataset(
            input_data,
            model_config,
            fasta_file,
            genotype_file,
        )

        attrs = model_wrapper.get_sequence_attributions(
            batch['ohe_seq'],
            batch['embed'],
            n_shuffles=n_shuffles,
            print_convergence_deltas=print_convergence_deltas
        ).squeeze(0).numpy()

        data.matrix = self.align_matrix_to_interval(
            matrix=attrs,
            matrix_interval=interval,
            target_interval=data.interval,
        ).T
        data.sequence_weights = data.matrix.sum(axis=1)
        return data
 
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
