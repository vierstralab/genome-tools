# genome-tools — notes for coding agents

Shared library of the Vierstra lab. Everything is keyed by `GenomicInterval`:
read data for an interval (`genome_tools.data.extractors`), manipulate intervals and
interval DataFrames (`genome_tools.genomic_interval`), plot an interval
(`genome_tools.plotting`, `genome_tools.plotting.modular_plot`).

Downstream packages (hotspot3, vinson, nmf-tools, footprint-tools) import it, so keep
public names and signatures stable, and keep the core light (see "Dependencies").

## Where things are

| Need | Look in |
|---|---|
| Interval class, DataFrame <-> intervals | `genome_tools/genomic_interval.py` |
| Read FASTA / BigWig / figwig / tabix / parquet / VCF for an interval | `genome_tools/data/extractors.py` |
| Open local/remote/gz files, write tabix | `genome_tools/helpers.py` |
| Backed AnnData from zarr | `genome_tools/data/anndata.py` |
| Single plot functions (`signal_plot`, `gene_annotation_plot`, `ideogram_plot`, `seq_plot`, `plot_motif_logo`) | `genome_tools/plotting/` |
| Axis helpers (`format_axes_to_interval`, `clear_spines`, `cm2inch`) | `genome_tools/plotting/utils.py` |
| Multi-track figures | `genome_tools/plotting/modular_plot/` (user guide: `README.md`, writing components: `README_develop.md`) |
| Signal smoothing, one-hot encoding | `genome_tools/utils/` |
| Runnable examples on public data | `examples/UsageExamples.ipynb` (data extraction, basic plots); `examples/README.md` indexes it and links component usage in vinson and regulome |

`genome_tools/data/dataset.py`, `loaders.py`, `samplers.py` are legacy training dataloaders
and will be removed. Don't use or extend them.

## GenomicInterval

- 0-based, half-open (BED). `to_ucsc()` prints the same numbers (`chr1:100-200`), no +1.
- Extra fields are kwargs: `GenomicInterval('chr1', 100, 200, name='peak1', score=5)`; they are carried by `copy`, `widen`, `shift`, `zoom`.
- Methods return a new interval unless `inplace=True`.
- `zoom(k)`: k > 0 zooms **in** (interval k times shorter), k < 0 zooms **out** (`zoom(-5)` is 5x longer).
- `resize(n)`: window of length n around the center. `center`: 1 bp interval at the middle.
- `widen` clips start at 0 and doesn't know chromosome sizes.
- `GenomicInterval.from_ucsc('chr1:1,000-2,000')` accepts commas.
- `VariantInterval(chrom, pos, pos + 1, ref, alt)` is a 1 bp interval; `widen`/`shift`/`zoom` on it raise, convert with `.gi()` first.
- `df_to_genomic_intervals(df, interval=None, extra_columns=())` expects `#chr`/`chrom`, `start`, `end` columns; passing `interval` keeps only rows overlapping it.

## Extractors

All extractors are context managers indexed by an interval:

```python
with FastaExtractor(fasta_file) as fasta:
    seq = fasta[interval].upper()           # str, case as in the file
with BigwigExtractor(bw_file) as bw:
    values = bw[interval]                   # np.ndarray, NaN where no data
with FigwigExtractor(fw_file, tracks='AG80080') as fw:
    values = fw[interval]                   # 1D for one track name, (n_tracks, len) otherwise
with TabixExtractor(bed_gz, columns=[...]) as tbx: # if columns not provided - infered from tabix header
    df = tbx[interval]                      # pd.DataFrame
```

- File handles (pysam, pyBigWig) are not thread/process safe: open one extractor per worker.
- `FigwigExtractor` imports `figwig` lazily (python >= 3.10, `pip install "genome-tools[figwig]"`).

## Modular plotting framework

- `IntervalPlotComponent` draws one row; it declares the data it needs via loaders
  (`@uses_loaders(...)` or `Component.with_loaders(..., new_class_name=..., doc=...)`).
- `PlotDataLoader._load(self, data, <args>)` reads files for `data.interval` and sets fields on
  `data` (a `DataBundle`). Its argument names become keyword arguments of the component.
- `IntervalPlotter(components, **loader_kwargs)` stacks components vertically.
  `plotter.plot(interval)` returns `(data, axes)`; or call `get_interval_data` then `plot_interval`.
  `interval` can be a dict `{interval_key: GenomicInterval}` to give rows different intervals.
- Loader argument priority (low -> high): loader defaults, `IntervalPlotter(**kwargs)`,
  component kwargs, `get_interval_data(**kwargs)`, loader-specific overrides on the component
  (`Component(SomeLoader={'arg': value})`), loader-specific overrides in `get_interval_data`.
- Discover a component's arguments: `help(Component)` or `inspect.signature(Component)`;
  required loader args show as `Required loader arg`. A missing one raises at load time naming
  the component and loader.
- Two components of the same class need distinct `name=`.
- `axes` (and `data`) are namedtuples keyed by component name: `axes.TrackComponent` or
  `getattr(axes, name)`, not `axes['TrackComponent']`.
- To plot a new data source with an existing plot, use `Component.with_loaders(..., doc=...)`
  instead of subclassing; subclass only to change `_plot`.
- To reuse a loader for another field, subclass it and override the attribute it sets
  (`SignalLoader.__signal_attr__`, `SegmentsLoader.__intervals_attr__`), e.g.
  `BackgroundBigwigLoader(SignalLoader)` sets `fit_threshold`. Give the subclass's `_load` its own
  argument name (`background_bw`, not `signal_file`): loader args are matched by name across all
  loaders of a component.
- hotspot3 components: `DensityAndBackgroundComponent` needs only the release bigWigs
  (`{id}.normalized_density.bw`, `{id}.background.bw`); `SignalAndMeanBGComponent`,
  `SmoothedSignalAndMeanBGComponent` and `FdrComponent` need parquet files from hotspot3's
  `debug/` output, which is not in the regulome release. 

## Dependencies

Core (`genome_tools`, `genomic_interval`, `helpers`, `data.extractors`, `plotting`,
`modular_plot.plot_components.basic`) must only need what is in `pyproject.toml`
dependencies. hotspot3 imports the core, so never import vinson, torch or other heavy
packages there, not even indirectly.

pip extras (`pyproject.toml`): `figwig` (`FigwigExtractor`), `anndata` (anndata, zarr 2, dask; only for
`genome_tools.data.anndata.read_zarr_backed`).

Optional dependencies are imported at module level by these plot component modules (import
the module only if the package is installed):

| Module (`modular_plot.plot_components.*`) | Needs |
|---|---|
| `basic` | — |
| `sequence` | bx-python |
| `variant` | bx-python, footprint-tools, hotspot3 |
| `footprint` | bx-python, footprint-tools, anndata |
| `hotspot3` | hotspot3 |
| `dhs` | nmf-tools |
| `prediction` | vinson (torch), anndata, bx-python |

## Working in this repo

- Python >= 3.8 for the core; match the style of the surrounding code (numpy-style or short
  docstrings, no type-annotation overhaul).
- Component docstrings follow the pattern in `plot_components/basic.py` and
  `plot_components/prediction.py`: what it plots, `Loaders:`, required/optional loader args,
  plot kwargs, a short example, `Returns:`. Keep the tables in `README.md` and
  `modular_plot/README.md` in sync when adding components.
- Docs and examples:
  - Examples use public files (regulome.altius.org, resources.altius.org), so they run anywhere;
  - `examples/UsageExamples.ipynb` covers data extraction and basic plot functions only;
    Don't repeat the same example in several places; link instead.

- A stale `build/` directory can shadow the source when building wheels; build from a clean
  checkout if packaging looks wrong.
