# Modular Plot Framework
The modular plot framework provides a flexible, composable system for creating complex genomic interval visualizations. It separates data loading from plotting logic, enabling reusable components that can be easily combined and configured. Details on individual components along with examples can be found in PlotComponents section

# Key concepts
- **IntervalPlotter**: The main orchestrator that manages components and creates figures

- **IntervalPlotComponent**: A plot component that renders data for a genomic interval

- **PlotDataLoader**: A data loader that fetches and processes data for components

- **DataBundle**: A container that holds data passed between loaders and components

# Basic functionality
## Quick start

```python
from genome_tools import GenomicInterval
from genome_tools.plotting.modular_plot import IntervalPlotter
from genome_tools.plotting.modular_plot.plot_components.basic import (
    IdeogramComponent,
    GencodeComponent,
    TrackComponent,
)
from genome_tools.plotting.ideogram import read_ideogram

# Define the genomic interval to visualize
interval = GenomicInterval('chr1', 1_000_000, 1_400_000)

# Define plot components with their heights and margins
plot_components = [
    IdeogramComponent(height=0.1, margins=(0.1, 0.1)),
    GencodeComponent(height=0.2, margins=(0.1, 0.1)),
    TrackComponent(height=1.5, margins=(0.1, 0.1)),
]

# Create the plotter
plotter = IntervalPlotter(plot_components, width=12)

# Plot in one step (load data + plot)
data, axes = plotter.plot(
    interval,
    loader_kwargs=dict(
        signal_file='path/to/signal.bw',  # Required by TrackComponent
        gencode_annotation_file='path/to/gencode.gtf',  # Required by GencodeComponent
        ideogram_data=read_ideogram('path/to/ideogram.txt'),  # Required by IdeogramComponent (pre-loaded object)
    )
)
```

## Two-step workflow (load then plot)

For more control, separate data loading from plotting:

```python
# Step 1: Load data
data = plotter.get_interval_data(
    interval,
    signal_file='path/to/signal.bw',
    gencode_annotation_file='path/to/gencode.gtf',
    ideogram_data=read_ideogram('path/to/ideogram.txt'),
)

# Step 2: Plot the data
axes = plotter.plot_interval(data, fig_width=12, height_scale=1.0)
```

This is useful when you want to either cache data for interactive exploration
or debug data loading separately from plotting

## Different intervals for components

Components can visualize different genomic intervals using `interval_key`:

```python
# Components with different interval keys
plot_components = [
    IdeogramComponent(height=0.1, interval_key='gene'),
    GencodeComponent(height=0.2, interval_key='gene'),
    TrackComponent(height=1.5, interval_key='dhs'),
]

plotter = IntervalPlotter(plot_components)

# Provide a dict of intervals
intervals = {
    'gene': GenomicInterval('chr1', 1_000_000, 1_400_000),
    'dhs': GenomicInterval('chr1', 1_150_000, 1_200_000),  # Zoomed in
}

data, axes = plotter.plot(
    intervals,
    signal_file='path/to/signal.bw',
    gencode_annotation_file='path/to/gencode.gtf',
    ideogram_data=read_ideogram('path/to/ideogram.txt'),
)
```

## Adding connectors

Connect components with lines or shaded areas to highlight specific regions:

```python
# Plot the interval
data, axes = plotter.plot(interval, **loader_kwargs)

# Add an area connector highlighting a regulatory region
plotter.plot_connector(
    axes,
    interval=GenomicInterval('chr1', 1_180_000, 1_190_000),
    start_component='GencodeComponent',
    end_component='TrackComponent',
    type='area',
    alpha=0.2,
    extend_to_top=False,
    extend_to_bottom=True,
)

# Add line connectors at specific positions
plotter.plot_connector(
    axes,
    positions=[1_185_000, 1_195_000, 1_205_000],
    start_component='GencodeComponent',
    end_component='TrackComponent',
    type='line',
    color='blue',
    linewidth=1,
)
```

**Connector types:**
- `'area'` - Shaded region between two x-coordinates (requires `interval` or `positions=[start, end]`)
- `'line'` - Vertical lines at specific positions (requires `positions=[x1, x2, ...]`, `interval` gets converted to `[x1, x2]`)

## (Not fully tested) Parallel Data Loading

For large datasets or slow loaders, enable parallel processing:

```python
data, axes = plotter.plot(
    interval,
    n_cpus=4,  # Use 4 processes for data loading
    **loader_kwargs
)
```

**Note:** Due to multiprocessing overhead, parallel loading may be slower for fast loaders. Best for I/O-bound operations.

# Available Plot Components

## Basic Components

### **TrackComponent**
Plots signal track from bigWig files.

**Loaders:** `SignalLoader`

**Required arguments:**
- `signal_file` (str): Path to bigwig file

**Example:**
```python
TrackComponent(
    height=1.5,
    margins=(0.1, 0.1),
    signal_file='accessibility.bw',
    color='blue',
    lw=1,
)
```

### **IdeogramComponent**
Displays chromosome ideogram with cytobands.

**Loaders:** `IdeogramLoader`

**Required arguments:**
- `ideogram_data` dict: Result of `genome_tools.plotting.ideogram.read_ideogram()`

**Example:**
```python
from genome_tools.plotting.modular_plot.plot_components.basic import IdeogramComponent
IdeogramComponent(
    height=0.1,
    margins=(0.1, 0.0),
    ideogram_data=read_ideogram('cytoBand.txt'),
)
```

### **GencodeComponent**
Plots gene annotations from GENCODE GTF files.

**Loaders:** `GencodeLoader`

**Required arguments:**
- `gencode_annotation_file` (str): Path to GENCODE GTF file

**Optional plot arguments:**
- `gene_symbol_exclude_regex` (str): Regex to exclude certain gene symbols (default: `r'^ENSG|^MIR|^LINC|.*-AS.*'`)

**Example:**
```python
from genome_tools.plotting.modular_plot.plot_components.basic import GencodeComponent
GencodeComponent(
    height=0.3,
    margins=(0.05, 0.05),
    gencode_annotation_file='gencode.v43.annotation.gtf',
    gene_symbol_exclude_regex=r'^ENSG',  # Exclude Ensembl IDs
)
```

## Sequence Components

### **SequencePlotComponent**
Plots DNA sequence as colored nucleotides.

**Loaders:** `FastaLoader`, `OHESequenceLoader`

**Required arguments:**
- `fasta_file` (str): Path to genome FASTA file

**Optional arguments:**
- `vocab` (dict): Dictionary for vocabulary ({letter: color}). Keys should go in order of ACGT. E.g. `vocab = {x: 'k' for x in ACGT}` to plot all letters as black

**Example:**
```python
from genome_tools.plotting.modular_plot.plot_components.sequence import SequencePlotComponent

SequencePlotComponent(
    height=0.2,
    fasta_file='hg38.fa',
)
```

### MotifHitsComponent
Displays transcription factor motif hits as sequence logos.

**Loaders:** `AnnotationRegionsLoader`, `MotifHitsLoader`, `MotifHitsSelectorLoader`

**Required arguments:**
- `annotation_regions` (List[GenomicInterval]): Regions to search for motifs
- `hits_file` (str): Path to motif hits tabix file. Has a `motif_id` column
- `motif_metadata` (pd.DataFrame): DataFrame with `motif_id` as index 

**Optional arguments:**
- `choose_by`: 
- `pack` (bool): Whether to pack overlapping motifs into rows

**Example:**
```python
from genome_tools.plotting.modular_plot.plot_components.sequence import MotifHitsComponent

MotifHitsComponent(
    height=0.5,
    annotation_regions=[GenomicInterval('chr1', 1_000_000, 1_001_000)],
    motif_annotations_path='motif_hits.bed.gz',  # tabix-indexed annotations
    motif_meta=motif_meta_df,  # DataFrame indexed by motif_id with metadata
    pack=True,
)
```

## DHS Components

### DHSIndexComponent
Displays DHS index regions as genomic segments.

**Loaders:** `DHSIndexLoader`

**Required arguments:**
- `dhs_index` (pd.DataFrame): DHS index rows within interval

**Example:**
```python
from genome_tools.plotting.modular_plot.plot_components.dhs import DHSIndexComponent

DHSIndexComponent(
    height=0.3,
    dhs_index=dhs_index_df,
)
```

### NMFTracksComponent (not used)
Plots multiple NMF component tracks stacked vertically.

**Loaders:** `ComponentTracksLoader`

**Required arguments:**
- `cutcounts_files` (dict[int, list[str]]): Map component index -> list of bigWig paths
- `nmf_components` (List[int]): Component indices to plot

**Optional plot arguments:**
- `common_lim` (bool): Use common y-axis limits across tracks
- `component_data` (pd.DataFrame): Metadata with 'index', 'color', 'name' (passed as plot kwarg)

**Example:**
```python
from genome_tools.plotting.modular_plot.plot_components.dhs import NMFTracksComponent

NMFTracksComponent(
    height=1.0,
    nmf_components=[0, 1, 2],
    cutcounts_files={0: files0, 1: files1, 2: files2},
    component_data=metadata_df,  # plot kwarg
)
```

# Using the same component multiple times

By default, components are named after their class. Multiple components within one plot can't have the same name. Customize names to use multiple instances:

```python
plotter = IntervalPlotter([
    TrackComponent(name='ATAC', signal_file='atac.bw', height=1.0),
    TrackComponent(name='H3K27ac', signal_file='h3k27ac.bw', height=1.0),
    TrackComponent(name='H3K4me3', signal_file='h3k4me3.bw', height=1.0),
])

data, axes = plotter.plot(interval)

# Access by custom name
axes.ATAC
axes.H3K27ac
axes.H3K4me3
```

# Debugging Data Loading

The `DataBundle` tracks which loaders have processed it:

```python
data = plotter.get_interval_data(interval, **loader_kwargs)

# Check loaded data for first component
component_data = data.IdeogramComponent
print(component_data)  # Shows all attached attributes
print(component_data.processed_loaders)  # Shows loader execution order

# Access specific attributes
print(component_data.interval)
print(component_data.signal)
```

# Adding additional objects on axes
Use returned axes object