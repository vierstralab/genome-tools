# Examples

- [UsageExamples.ipynb](UsageExamples.ipynb): genomic intervals, reading FASTA / BigWig / tabix data for an interval, and the basic plot functions (ideogram, motif logo, signal, gene annotation). Uses public data only.

To plot several tracks together (`IntervalPlotter` and plot components, including vinson predictions and hotspot3 results), see the [modular plotting guide](../genome_tools/plotting/modular_plot/README.md).

Additional examples that use modular plot framework:
- [hotspot3 figures](https://github.com/vierstralab/regulome/blob/master/Figure1.index/hotspot3.ipynb) in regulome: cut counts with the fitted background, peak delineation

- [vinson tutorials](https://github.com/vierstralab/vinson/tree/main/tutorials): [predicted vs true signal tracks](https://github.com/vierstralab/vinson/blob/main/tutorials/2.Predict_density.ipynb), [attribution logos and motif hits](https://github.com/vierstralab/vinson/blob/main/tutorials/3.Interpretation.ipynb), [scoring a variant](https://github.com/vierstralab/vinson/blob/main/tutorials/4.Predict_variant.ipynb) (`FastaExtractor`, `read_zarr_backed`)
- regulome figure notebooks with the vinson model:
  - [per-sample DHS validation](https://github.com/vierstralab/regulome/blob/master/Figure7.SequenceModel/SeqModel.PerSampleDHSValidation.ipynb): predicted vs true signal with peaks
  - [cell-selective DHS validation](https://github.com/vierstralab/regulome/blob/master/Figure7.SequenceModel/SeqModel.CellSelectiveValidation.ipynb): predicted signal, attributions, motif hits and footprints
  - [between-species predictions](https://github.com/vierstralab/regulome/blob/master/Figure8.SpeciesPrediction/SeqModel.BetweenSpeciesPredictions.ipynb): attributions and motif hits aligned between species
  - [mouse cell-selective tracks](https://github.com/vierstralab/regulome/blob/master/Figure8.SpeciesPrediction/SeqModel.CellSelective.tracks_mouse.ipynb): predicted signal for mouse samples
