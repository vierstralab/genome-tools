# genome-tools

genome-tools is a shared Python toolkit built around a single core object, the `GenomicInterval`. Everything else is keyed by it:

- **Data extraction and manipulation** — pull sequence, signal and records for an interval from FASTA, BigWig, FigWig, tabix, parquet, D4 and VCF/BAM files, and convert between intervals and pandas DataFrames.
- **Plotting** — standalone plot functions and a modular framework for composing multi-track, genome-browser-style figures of an interval.

It underlies several downstream Vierstra lab packages as well (e.g. hotspot3, vinson, nmf-tools).

Examples: [UsageExamples.ipynb](examples/UsageExamples.ipynb) (reading data and basic plots, runs on public data); more, including multi-track plots and their use in vinson and regulome, in [examples/README.md](examples/README.md).

Coding agents: read [AGENTS.md](AGENTS.md) first.

## Table of contents

- [Repository layout](#repository-layout)
- [Installation](#installation)
- [Core components](#core-components)
  - [Genomic intervals](#genomic-intervals)
  - [Overlaps between interval tables](#overlaps-between-interval-tables)
  - [Data extraction](#data-extraction)
  - [Plotting](#plotting)
  - [Modular plotting framework](#modular-plotting-framework)
- [Authors](#authors)

## Repository layout

- [`genome_tools/genomic_interval.py`](genome_tools/genomic_interval.py) — `GenomicInterval` / `VariantInterval` classes and DataFrame <-> interval conversion helpers.
- [`genome_tools/bed_overlaps.py`](genome_tools/bed_overlaps.py) — overlap flags and distance to the closest interval between BED-like DataFrames (`pyranges` extra).
- [`genome_tools/data/`](genome_tools/data/) — file extractors (FASTA, BigWig, FigWig, tabix, parquet, D4, VCF genotypes/reads) and `AnnData` helpers.
- [`genome_tools/plotting/`](genome_tools/plotting/) — standalone plot functions (`signal_plot`, `gene_annotation_plot`, `ideogram_plot`, sequence/motif-logo plotting, color palettes/colormaps) and the [modular plotting framework](genome_tools/plotting/modular_plot/) for building composable, multi-track interval plots.
- [`genome_tools/utils/`](genome_tools/utils/) — signal smoothing/aggregation and sequence helpers.
- [`examples/`](examples/) — [UsageExamples.ipynb](examples/UsageExamples.ipynb) and an [index of examples](examples/README.md), including plot components used in vinson and regulome.


## Installation

### Prerequisites
`python >= 3.8`

### Install genome-tools
```bash
pip install genome-tools
```

```bash
git clone https://github.com/vierstralab/genome-tools.git
cd genome-tools
pip install .
```
Or directly from git:
```bash
pip install git+https://github.com/vierstralab/genome-tools.git
```

Optional extras:
```bash
pip install "genome-tools[figwig]"   # FigwigExtractor / FigwigTrackComponent, python >= 3.10
pip install "genome-tools[anndata]"  # read_zarr_backed (anndata, zarr 2, dask)
pip install "genome-tools[sampling]" # genome_tools.utils.sampling (numba)
pip install "genome-tools[pyranges]" # genome_tools.bed_overlaps
```
Some plot components need the package they visualize (`vinson` for prediction/attribution components, `hotspot3`, `footprint-tools`, `bx-python` for between-species alignment). They are imported only when you import those component modules.

## Core components

### Genomic intervals

`GenomicInterval` ([`genomic_interval.py`](genome_tools/genomic_interval.py)) is a BED-style coordinate object that supports arbitrary extra fields (e.g. `score`, `name`, custom metadata columns) attached as kwargs:

```python
from genome_tools import GenomicInterval

iv = GenomicInterval("chr1", 1_000_000, 1_001_000, name="peak1", score=10)
iv.to_ucsc()          # "chr1:1000000-1001000"
GenomicInterval.from_ucsc("chr1:1,000,000-1,001,000")

iv.widen(500)          # expand both ends by 500bp
iv.zoom(2.0)           # zoom in 2x around the center (zoom(-2.0) zooms out)
iv.shift(100)          # shift coordinates
iv.overlaps(other_iv)  # bool
```

`VariantInterval` is a 1bp specialization carrying `ref`/`alt`/`pos`; Helpers `df_to_genomic_intervals`, `df_to_variant_intervals`, `intervals_to_df`, and `filter_df_to_interval` convert between pandas DataFrames and lists of intervals.

### Overlaps between interval tables

[`genome_tools/bed_overlaps.py`](genome_tools/bed_overlaps.py) annotates a BED-like DataFrame (`#chr`/`chrom`, `start`, `end`) against other ones. Results are aligned to the index of the first DataFrame. Needs the `pyranges` extra.

```python
from genome_tools.bed_overlaps import overlaps_any, distance_to_closest

overlaps_any(peaks, tss)                   # bool Series: overlaps any TSS
overlaps_any(peaks, tss=tss, enh=enh)      # DataFrame, one bool column per map
distance_to_closest(peaks, tss)            # bp to the closest TSS; 0 if overlapping, NaN if none on the chromosome
```

### Data extraction

[`genome_tools/data/extractors.py`](genome_tools/data/extractors.py) provides context-manager extractors keyed by `GenomicInterval`, each wrapping the relevant file format:
- `FastaExtractor` / `InMemoryFastaExtractor` — pysam-backed and fully-in-RAM FASTA sequence access
- `BigwigExtractor` — per-base signal from BigWig
- `FigwigExtractor` — per-base signal from one or many tracks of a [figwig](https://github.com/jmschrei/figwig) (`.fw`) file, returns `(n_tracks, len(interval))` (or 1D for a single track name)
- `TabixExtractor` — tabix-indexed BED/GFF-style records
- `ChromParquetExtractor` — per-chromosome-partitioned parquet signal
- `D4Extractor` — D4 coverage files
- `VariantGenotypeExtractor` / `AllelicReadsExtractor` — genotypes and allele-specific read counts from VCF/BAM

⚠️ Most extractors wrap non-thread-safe file handles (`pysam`, `pyBigWig`) — open separate instances per worker when using multiprocessing.

[`genome_tools/data/anndata.py`](genome_tools/data/anndata.py) provides `read_zarr_backed()` for opening large sample-by-peak AnnData saved with `AnnData.write_zarr` without reading every matrix into RAM (older anndata versions can't do this). Needs the `anndata` extra (`pip install "genome-tools[anndata]"`). Latest versions of anndata already implement analogous `read_lazy` function

```python
from genome_tools.data.anndata import read_zarr_backed

adata = read_zarr_backed('index.anndata.zarr')               # lazy_attrs=('layers', 'varm')
adata.layers['counts']                                       # dask array, nothing read yet
adata[:, mask].layers['counts'].compute()                    # read only the selected columns
```

- `X`, `obs`, `var`, `obsm`, `obsp`, `varp` and `uns` are read into memory.
- Every entry of the groups named in `lazy_attrs` (default `layers` and `varm`) becomes a `dask.array`, read only when you call `.compute()`. A sparse entry in one of these groups is read into memory instead.
- Tested with anndata 0.10.9, zarr 2.18

### Plotting

[`genome_tools/plotting/`](genome_tools/plotting/__init__.py) exposes standalone plotting functions: `signal_plot` (continuous tracks), `gene_annotation_plot` (Gencode-style gene models), `ideogram_plot` (chromosome ideograms), sequence/motif-logo plotting, plus axis helpers (`format_axes_to_interval`, `pack_rows`, `clear_spines`) and color palettes/colormaps.

### Modular plotting framework

[`genome_tools/plotting/modular_plot/`](genome_tools/plotting/modular_plot/) decouples **loading** data for a genomic interval from **rendering** it (`DataBundle`, `PlotDataLoader`, `IntervalPlotComponent`, `IntervalPlotter`), so multi-track figures (ideogram + gene model + signal + motifs + predictions, etc.) can be assembled from reusable, mix-and-match components instead of per-figure plotting code.

```python
from genome_tools import GenomicInterval
from genome_tools.plotting.ideogram import read_ideogram
from genome_tools.plotting.modular_plot import IntervalPlotter
from genome_tools.plotting.modular_plot.plot_components.basic import (
    IdeogramComponent, GencodeComponent, TrackComponent,
)

plotter = IntervalPlotter(
    [
        IdeogramComponent(height=0.1, margins=0),
        GencodeComponent(height=0.3),
        TrackComponent(signal_file="https://regulome.altius.org/data/per_sample/AG95384/AG95384.normalized_density.bw", height=0.5, color="tab:blue"),
    ],
    ideogram_data=read_ideogram("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz"),   # loader args shared by components
    gencode_annotation_file="https://resources.altius.org/~sabramov/collaborations/vinson/examples/gencode.v49.basic.annotation.sorted.gff3.gz",
)
data, axes = plotter.plot(GenomicInterval.from_ucsc("chr14:72,600,000-72,800,000"))
```

Commonly used components:

| Component | Module | Shows | Key loader args |
|---|---|---|---|
| `IdeogramComponent` | `plot_components.basic` | chromosome ideogram | `ideogram_data` |
| `GencodeComponent` | `plot_components.basic` | gene models | `gencode_annotation_file` |
| `TrackComponent` | `plot_components.basic` | BigWig signal | `signal_file` |
| `SequencePlotComponent` | `plot_components.sequence` | DNA sequence | `fasta_file` |
| `PredictedSignalComponent` | `plot_components.prediction` | vinson predicted signal | `model_wrapper`, `sample_id`, `embeddings`, `fasta_file` |
| `AttributionsFromBatchComponent` | `plot_components.prediction` | vinson attributions as a sequence logo | `model_wrapper`, `batch` |
| `AttributionsWeightedMotifHitsFromBatchComponent` | `plot_components.prediction` | motif logos in annotated regions | + `annotation_regions`, `motif_annotations_path`, `motif_meta` |

Every component lists its loader arguments (required and with defaults) in its signature, and describes them in its docstring: `help(TrackComponent)` or `?TrackComponent` in jupyter.

See [`genome_tools/plotting/modular_plot/README.md`](genome_tools/plotting/modular_plot/README.md) for the full guide — quick start, available plot components, overriding loader parameters, and debugging data loading — and [`README_develop.md`](genome_tools/plotting/modular_plot/README_develop.md) for writing your own components.


## Authors

Developed by (alphabetically by last name):

- Sergey Abramov
- Alexandr Boytsov
- Jeff Vierstra

Altius Institute for Biomedical Sciences
