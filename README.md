# genome-tools

genome-tools is a shared Python toolkit for working with genomic intervals and sequencing data: interval arithmetic, BED/tabix/BigWig/FASTA I/O, batched data loading for model training, and a modular framework for composing multi-track genome-browser-style plots. It underlies several downstream Altius/Vierstra-lab packages (e.g. hotspot3, vinson, nmf-tools) rather than being a standalone end-user tool.

## Table of contents

- [Repository layout](#repository-layout)
- [Installation](#installation)
- [Core components](#core-components)
  - [Genomic intervals](#genomic-intervals)
  - [Data extraction and loading](#data-extraction-and-loading)
  - [Plotting](#plotting)
  - [Modular plotting framework](#modular-plotting-framework)
- [Authors](#authors)

## Repository layout

- [`genome_tools/genomic_interval.py`](genome_tools/genomic_interval.py) — `GenomicInterval` / `VariantInterval` classes and DataFrame <-> interval conversion helpers.
- [`genome_tools/data/`](genome_tools/data/) — file extractors (FASTA, BigWig, tabix, parquet, D4, VCF genotypes/reads),
- [`genome_tools/plotting/`](genome_tools/plotting/) — standalone plot functions (`signal_plot`, `gene_annotation_plot`, `ideogram_plot`, sequence/motif-logo plotting, color palettes/colormaps) and the [modular plotting framework](genome_tools/plotting/modular_plot/) for building composable, multi-track interval plots.
- [`genome_tools/utils/`](genome_tools/utils/) — signal smoothing/aggregation and sequence helpers.


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

## Core components

### Genomic intervals

`GenomicInterval` ([`genomic_interval.py`](genome_tools/genomic_interval.py)) is a BED-style coordinate object that supports arbitrary extra fields (e.g. `score`, `name`, custom metadata columns) attached as kwargs:

```python
from genome_tools import GenomicInterval

iv = GenomicInterval("chr1", 1_000_000, 1_001_000, name="peak1", score=10)
iv.to_ucsc()          # "chr1:1000000-1001000"
GenomicInterval.from_ucsc("chr1:1,000,000-1,001,000")

iv.widen(500)          # expand both ends by 500bp
iv.zoom(2.0)           # zoom out 2x around the center
iv.shift(100)          # shift coordinates
iv.overlaps(other_iv)  # bool
```

`VariantInterval` is a 1bp specialization carrying `ref`/`alt`/`pos`; Helpers `df_to_genomic_intervals`, `df_to_variant_intervals`, `intervals_to_df`, and `filter_df_to_interval` convert between pandas DataFrames and lists of intervals.

### Data extraction and loading

[`genome_tools/data/extractors.py`](genome_tools/data/extractors.py) provides context-manager extractors keyed by `GenomicInterval`, each wrapping the relevant file format:
- `FastaExtractor` / `InMemoryFastaExtractor` — pysam-backed and fully-in-RAM FASTA sequence access
- `BigwigExtractor` — per-base signal from BigWig
- `TabixExtractor` — tabix-indexed BED/GFF-style records
- `ChromParquetExtractor` — per-chromosome-partitioned parquet signal
- `D4Extractor` — D4 coverage files
- `VariantGenotypeExtractor` / `AllelicReadsExtractor` — genotypes and allele-specific read counts from VCF/BAM

⚠️ Most extractors wrap non-thread-safe file handles (`pysam`, `pyBigWig`) — open separate instances per worker when using multiprocessing.

[`genome_tools/data/anndata.py`](genome_tools/data/anndata.py) provides `read_zarr_backed()` for loading large sample-by-peak matrices as `AnnData` objects for legacy versions of anndata package.

### Plotting

[`genome_tools/plotting/`](genome_tools/plotting/__init__.py) exposes standalone plotting functions: `signal_plot` (continuous tracks), `gene_annotation_plot` (Gencode-style gene models), `ideogram_plot` (chromosome ideograms), sequence/motif-logo plotting, plus axis helpers (`format_axes_to_interval`, `pack_rows`, `clear_spines`) and color palettes/colormaps.

### Modular plotting framework

[`genome_tools/plotting/modular_plot/`](genome_tools/plotting/modular_plot/) decouples **loading** data for a genomic interval from **rendering** it (`DataBundle`, `PlotDataLoader`, `IntervalPlotComponent`, `IntervalPlotter`), so multi-track figures (ideogram + gene model + signal + motifs + predictions, etc.) can be assembled from reusable, mix-and-match components instead of per-figure plotting code.

See [`genome_tools/plotting/modular_plot/README.md`](genome_tools/plotting/modular_plot/README.md) for the full guide — quick start, available plot components, overriding loader parameters, and debugging data loading — and [`README_develop.md`](genome_tools/plotting/modular_plot/README_develop.md) for writing your own components.


## Authors

Developed by (alphabetically by last name):

- Sergey Abramov
- Alexandr Boytsov
- Jeff Vierstra

Altius Institute for Biomedical Sciences
