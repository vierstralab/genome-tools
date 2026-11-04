import numpy as np
import anndata as ad
import pandas as pd

from scipy.interpolate import interp1d

from typing import Tuple

from vinson.inference import assemble_batch, assemble_batch_from_data
from vinson.utils.data_formatting import extract_data_from_backed_anndata, VinsonData
from vinson.postprocessing.interpretation import ModelWrapper

from genome_tools import GenomicInterval
from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle
from genome_tools.data.utils import realign_matrix
from genome_tools.plotting.modular_plot.helpers.species_mapper import map_matrix_to_interval, BetweenSpeciesMap


class BatchLoader(PlotDataLoader):
    """Use a pre-built model input `batch` as is. Sets `data.batch`."""
    required_loader_kwargs=['batch']


class BatchFromIntervalCenterLoader(PlotDataLoader):
    """Build a single-example vinson batch centered on the interval center for `sample_id`. Sets `data.batch`."""

    def _load(self, data: DataBundle,
                    sample_id: str,
                    embeddings: pd.DataFrame,
                    fasta_file: str,
                    genotype_file: str=None,
            ):

        batch = assemble_batch(
            sample_id,
            data.interval,
            embeddings,
            fasta_file,
            genotype_file=genotype_file,
        )

        data.batch = batch
        return data


class BatchFromSteppedIntervalLoader(PlotDataLoader):
    """Build a vinson batch with model inputs centered every `step` bp across the interval for `sample_id`. Sets `data.batch`."""
    def _load(self, data: DataBundle,
                sample_id: str,
                embeddings: pd.DataFrame,
                fasta_file: str,
                genotype_file: str=None,
                step=20
        ):
        batch = assemble_batch(
            sample_id,
            data.interval,
            embeddings,
            fasta_file,
            genotype_file=genotype_file,
            step=step,
        )

        data.batch = batch
        return data


class PredictedSignalLoader(PlotDataLoader):
    """Run `model_wrapper` on `data.batch` and interpolate predictions to per-base `data.signal` over the interval."""
    def _load(self, data: DataBundle, model_wrapper: ModelWrapper, interp1d_kind='linear'):
        initial_interval: GenomicInterval = data.interval

        batch = data.batch 

        # replace with predict_with_model

        predictions_coords = batch['summit']
        X_seq = batch["ohe_seq"].to(model_wrapper.model.device, non_blocking=True)
        X_embed = batch["embed"].to(model_wrapper.model.device, non_blocking=True)

        pred_density = model_wrapper(X_seq, X_embed).detach().cpu().numpy().squeeze()

        full_positions = np.arange(initial_interval.start, initial_interval.end)

        order = np.argsort(predictions_coords)

        data.signal = interp1d(
            np.asarray(predictions_coords)[order],
            pred_density[order],
            kind=interp1d_kind,
            assume_sorted=True,
        )(full_positions)
        return data


class AttributionsLoader(PlotDataLoader):
    """Compute DeepLIFT/SHAP attributions for a single-example `data.batch`.

    Sets `data.raw_attributions` (4 x model input length) and `data.raw_attributions_coords` (model input interval).
    """
    def _load(self, data: DataBundle,
                model_wrapper: ModelWrapper,
                print_convergence_deltas=False,
                n_shuffles=20,
                random_state=42,
        ):
        batch = data.batch
        target_interval: GenomicInterval = data.interval

        assert len(batch['chrom']) == 1
        batch_interval = GenomicInterval(
            batch['chrom'][0],
            batch['summit'][0],
            batch['summit'][0] + 1
        ).resize(
            batch['ohe_seq'][0].shape[1]
        )
        assert target_interval.overlaps(batch_interval), 'Interval in loaded batch does not overlap target interval. Check batch coordinates'

        attrs = model_wrapper.get_sequence_attributions(
            batch['ohe_seq'],
            batch['embed'],
            n_shuffles=n_shuffles,
            print_convergence_deltas=print_convergence_deltas,
            random_state=random_state
        ).squeeze(0).numpy()
        data.raw_attributions = attrs
        data.raw_attributions_coords = batch_interval
        return data


class AlignedAttributionsLoader(PlotDataLoader):
    """Crop/pad raw attributions to the interval. Sets `data.matrix` and `data.sequence_weights`."""
    def _load(self, data: DataBundle):
        data.matrix = realign_matrix(
            matrix=data.raw_attributions,
            matrix_interval=data.raw_attributions_coords,
            target_interval=data.interval,
        ).T
        data.sequence_weights = data.matrix.sum(axis=1)
        return data

 
class BetweenSpeciesAlignedAttributionsLoader(PlotDataLoader):
    """Map raw attributions onto the interval through a `BetweenSpeciesMap`. Sets `data.matrix` and `data.sequence_weights`."""
    def _load(self, data: DataBundle, mapping: BetweenSpeciesMap):
        data.matrix = map_matrix_to_interval(
            matrix=data.raw_attributions,
            matrix_interval=data.raw_attributions_coords,
            target_interval=data.interval,
            mapping=mapping
        ).T
        data.sequence_weights = data.matrix.sum(axis=1)
        return data


# Load batch from anndata, defunc
class BatchFromAnndataLoader(PlotDataLoader):

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

    def _load(self, data: DataBundle,
                sample_id: str,
                anndata: ad.AnnData,
                dhs_id: str,
                fasta_file: str,
                genotype_file: str=None,
        ):
        input_data, interval = self.from_backed_anndata(
            anndata,
            sample_id=sample_id,
            dhs_id=dhs_id,
        )
        assert interval.overlaps(data.interval), f"Data interval {data.interval} does not overlap with the dhs interval {interval}, {dhs_id}"

        batch = assemble_batch_from_data(
            input_data,
            fasta_file,
            genotype_file=genotype_file,
        )
        data.batch = batch
        return data
