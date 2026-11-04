# Examples

- [UsageExamples.ipynb](UsageExamples.ipynb): genomic intervals, reading FASTA / BigWig / tabix data for an interval, single plots (ideogram, motif logo, signal, gene annotation), multi-track figures with `IntervalPlotter`, zoom levels with connectors, and writing your own plot component (a heatmap). Uses public data only.

Examples in the [vinson](https://github.com/vierstralab/vinson) tutorials (need `vinson` and a trained model to run):

- [2.Predict_density](https://github.com/vierstralab/vinson/blob/main/tutorials/2.Predict_density.ipynb): true vs predicted DNase I signal tracks under an ideogram and gene annotation (`PredictedSignalComponent`).
- [3.Interpretation](https://github.com/vierstralab/vinson/blob/main/tutorials/3.Interpretation.ipynb): attribution logos, motif hits in annotated regions and connectors (`AttributionsFromBatchComponent`, `AttributionsWeightedMotifHitsFromBatchComponent`).
- [4.Predict_variant](https://github.com/vierstralab/vinson/blob/main/tutorials/4.Predict_variant.ipynb): `FastaExtractor` and `read_zarr_backed` to score a variant with the model.
