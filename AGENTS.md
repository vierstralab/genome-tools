# genome-tools — notes for coding agents

Shared library of the Vierstra lab. Everything is keyed by `GenomicInterval`:
read data for an interval (`genome_tools.data.extractors`), manipulate intervals and
interval DataFrames (`genome_tools.genomic_interval`), plot an interval
(`genome_tools.plotting`, `genome_tools.plotting.modular_plot`).

When to use: code that imports `genome_tools`, or works with genomic intervals: building or
shifting/zooming intervals, reading FASTA / bigWig / figwig / tabix data for a region,
overlapping BED-like DataFrames, plotting a region (signal tracks, genes, ideogram, motif logos),
or making multi-track figures with `IntervalPlotter` and plot components (also from
vinson, hotspot3 or regulome notebooks).

Downstream packages (hotspot3, vinson, nmf-tools, footprint-tools) import it, so keep
public names and signatures stable, and keep the core light (see "Dependencies").

## Where things are

| Need | Look in |
|---|---|
| Interval class, DataFrame <-> intervals | `genome_tools/genomic_interval.py` |
| Read FASTA / BigWig / figwig / tabix / parquet / VCF for an interval | `genome_tools/data/extractors.py` |
| Open local/remote/gz files, write tabix | `genome_tools/helpers.py` |
| Overlap flags / distance to closest between BED-like DataFrames (`overlaps_any`, `distance_to_closest`, `to_pyranges`) | `genome_tools/bed_overlaps.py` |
| Backed AnnData from zarr | `genome_tools/data/anndata.py` |
| Single plot functions (`signal_plot`, `gene_annotation_plot`, `ideogram_plot`, `seq_plot`, `plot_motif_logo`) | `genome_tools/plotting/` |
| Axis helpers (`format_axes_to_interval`, `clear_spines`, `cm2inch`) | `genome_tools/plotting/utils.py` |
| Multi-track figures | `genome_tools/plotting/modular_plot/` (user guide: `README.md`, writing components: `README_develop.md`, notebooks: `examples/`) |
| Signal smoothing, one-hot encoding | `genome_tools/utils/` |
| Runnable examples on public data | `examples/UsageExamples.ipynb` (data extraction, basic plots); `examples/README.md` indexes it and links component usage in vinson and regulome |

`genome_tools/data/footprints/`: legacy dataloaders used by footprint-tools.

## GenomicInterval

- 0-based, half-open (BED). `to_ucsc()` prints UCSC formatted string (`chr1:100-200`)
  `to_ucsc(commas=True)` adds thousands separators (`chr1:1,000-2,000`).
- Extra fields are kwargs: `GenomicInterval('chr1', 100, 200, name='peak1', score=5)`; they are carried by `copy`, `widen`, `shift`, `zoom`.
- Methods return a new interval unless `inplace=True`.
- `zoom(k)`: k > 0 zooms **in** (interval k times shorter), k < 0 zooms **out** (`zoom(-5)` is 5x longer).
- `resize(n)`: window of length n around the center. `center`: 1 bp interval at the middle.
- `widen` clips start at 0 and doesn't know chromosome sizes.
- `GenomicInterval.from_ucsc('chr1:1,000-2,000')` accepts commas.
- `VariantInterval(chrom, pos, pos + 1, ref, alt)` is a 1 bp interval; `widen`/`shift`/`zoom` on it raise, convert with `.gi()` first.
- `df_to_genomic_intervals(df, interval=None, extra_columns=())` expects `#chr`/`chrom`, `start`, `end` columns; passing `interval` keeps only rows overlapping it.
  Only these columns and `extra_columns` are carried over.
- `df_to_variant_intervals(df, interval=None, extra_columns=())` also needs `ref` and `alt`; with a
  `VariantInterval` as `interval`, rows must match its `ref`/`alt` too (pass `.gi()` to filter by position only).
- `intervals_to_df(intervals)`, the inverse of `df_to_*_intervals`: DataFrame with `chrom`, `start`,
  `end` (+ `ref`, `alt`), then a column per extra field (always includes `name`). All intervals must
  be of the same class and have the same extra field names (in any order, values differ), as
  intervals made by `df_to_*_intervals` do.

## Interval tables (bed_overlaps)

- Input DataFrames need `#chr` or `chrom`, plus `start` and `end` (0-based, half-open); other columns are ignored.
- Both functions return one value per row of `ref_df`, with the same index as `ref_df`.
- `distance_to_closest` follows the pyranges/bedtools convention: intervals that touch end-to-start are 1 bp apart; the result is NaN if the map has no intervals on that chromosome.
- Use these on DataFrames; use `GenomicInterval.overlaps` on single intervals.

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
- Loaders must never edit an object in place (`data.signal[:] = 0`, `df['x'] = ...`,
  `inplace=True`), neither existing `data` fields nor loader args: assign a new object
  (`data.signal = data.signal * 0`). Even that should be rare: prefer setting a new field
  (`data.masked_signal = ...`) over replacing one another loader set. Then every field keeps the
  meaning given by the loader that set it, and loaders can be combined in any component
  (`with_loaders`) without depending on their order; `DataBundle` warns when a field is
  overridden. Components share these objects: `get_interval_data`
  (`use_loaders_cache=True`, default) runs a loader once for components with the same interval,
  previous loaders and loader args (same objects, compared by `id`), and loader args are passed
  as is. In-place edits can't be detected and silently change other components' data.
- `IntervalPlotter(components, **loader_kwargs)` stacks components vertically.
  `plotter.plot(interval)` returns `(data, axes)`; or call `get_interval_data` then `plot_interval`.
  `interval` can be a dict `{interval_key: GenomicInterval}` to give rows different intervals.
- `plot_interval(data, subplot_spec=outer_gs[0, 1])` draws into a cell of an existing figure's
  gridspec instead of a new figure (`fig_width`/`fig_height`/`height_scale` are then ignored).
  `plot()` always makes a new figure.
- `plotter.get_loading_workflow(interval, **loader_kwargs)` shows the loader calls
  `get_interval_data` would run (shared calls once) without loading anything; `.plot()` draws it
  as a tree, `.tree` is the `nx.DiGraph`.
- Loader argument priority (low -> high), init before runtime and the more specific winning within each:
  loader defaults, `IntervalPlotter(**kwargs)`, component kwargs, loader-specific overrides on the
  component (`Component(SomeLoader={'arg': value})`), `get_interval_data(**kwargs)`,
  `get_interval_data(ComponentName={...})`, loader-specific overrides in `get_interval_data`
  (`SomeLoader={...}`, or `ComponentName={'SomeLoader': {...}}`).
- In a notebook, a plotter, component, loaded data (namedtuple of `DataBundle`) or loading workflow
  shown as the last line of a cell (or `display(...)`) renders as markdown (`_repr_markdown_`):
  resolved loader args and where they were set, bundle fields, shared loader calls.
  `help()` / `?` show the signature and docstring, not this summary. Use it to debug a plot
  (what is set, what is missing, what was loaded); `modular_plot/examples/Debugging.ipynb` shows it all.
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
`modular_plot.plot_components.basic`): hotspot3 imports the core, so never import vinson,
torch or other heavy packages there, not even indirectly.

pip extras (`pyproject.toml`): `figwig` (`FigwigExtractor`), `anndata` (anndata, zarr 2, dask; only for
`genome_tools.data.anndata.read_zarr_backed`), `pyranges` (imported inside
`genome_tools.bed_overlaps`), `sampling` (numba, imported at module level by `genome_tools.utils.sampling`).

Optional dependencies of the plot component modules. "Needs at import" is imported at module
level (directly or through the `loaders/` and `helpers/` modules the component module imports):
import the module only if these are installed. "Needs when used" is imported inside a function,
only when that component/loader runs:

| Module (`modular_plot.plot_components.*`) | Needs at import | Needs when used |
|---|---|---|
| `basic` | — | figwig (`FigwigTrackComponent`) |
| `sequence` | bx-python | — |
| `variant` | bx-python | hotspot3 (`CutcountsComponent`, `AllelicCutcountsComponent`), footprint-tools (`ReadsLoader` in `CutcountsComponent`) |
| `footprint` | bx-python, footprint-tools, anndata | — |
| `hotspot3` | hotspot3 | — |
| `dhs` | — | nmf-tools (`DHSLoadingsComponent`) |
| `prediction` | vinson (torch), anndata, bx-python | — |

## Working in this repo

- Python >= 3.9 for the core; match the style of the surrounding code (numpy-style or short
  docstrings, no type-annotation overhaul).
- Component docstrings follow the pattern in `plot_components/basic.py` and
  `plot_components/prediction.py`: what it plots, `Loaders:`, required/optional loader args,
  plot kwargs, a short example, `Returns:`. Keep the tables in `README.md` and
  `modular_plot/README.md` in sync when adding components.
- Docs and examples:
  - Examples use public files (regulome.altius.org, resources.altius.org), so they run anywhere;
  - To check a change on real data (instead of mock arrays), read these files remotely; extractors
    and loaders accept the URLs directly, no download needed:
    - regulome per sample (`{s}` e.g. `AG95384`, `AG95451`), base `https://regulome.altius.org/data/per_sample/{s}/`:
      `{s}.normalized_density.bw`, `peaks/{s}.background.bw`,
      `peaks/fdr0.05/{s}.hotspots.fdr0.05.bed.gz`, `peaks/fdr0.05/{s}.peaks.fdr0.05.bed.gz`
    - GENCODE: `https://resources.altius.org/~sabramov/collaborations/vinson/examples/gencode.v49.basic.annotation.sorted.gff3.gz`
    - FASTA: `https://resources.altius.org/~sabramov/collaborations/vinson/examples/GRCh38_no_alts.fa`
    - ideogram: `read_ideogram('http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBandIdeo.txt.gz')`
    - region with signal for all of the above: `chr2:60450000-60540000` (BCL11A)
  - `examples/UsageExamples.ipynb` covers data extraction and basic plot functions only;
    modular plot examples are in `genome_tools/plotting/modular_plot/examples/`: `Debugging.ipynb`
    (notebook displays), `SetLoaderAndPlotKwargs.ipynb` (an example for each priority level of loader args and plot kwargs).
    Don't repeat the same example in several places; link instead.

- A stale `build/` directory can shadow the source when building wheels; build from a clean
  checkout if packaging looks wrong.
