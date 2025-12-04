# Modular Plot Framework

The modular plot framework provides a flexible, composable system for creating complex genomic interval visualizations. It separates data loading from plotting logic, enabling reusable components that can be easily combined and configured.

## Key Concepts

- **IntervalPlotter**: The main orchestrator that manages components and creates figures

- **IntervalPlotComponent**: A plot component that renders data for a genomic interval

- **PlotDataLoader**: A data loader that fetches and processes data for components

- **DataBundle**: A container that holds data passed between loaders and components

### Architecture

The framework follows a clean separation of concerns:

1. **Data Loaders** - Fetch and preprocess data (e.g., from bigWig, tabix, FASTA files)
2. **Plot Components** - Visualize the loaded data on matplotlib axes
3. **Interval Plotter** - Coordinates data loading and plotting across multiple components

This design allows you to:
- Mix and match components without modifying their code
- Share loaders across different components
- Override loader parameters at multiple levels (global, per-component, runtime)
- Create complex multi-panel figures with minimal boilerplate

## Basic Usage with IntervalPlotter

### Quick Start

```python
from genome_tools import GenomicInterval
from genome_tools.plotting.modular_plot import IntervalPlotter
from genome_tools.plotting.modular_plot.plot_components.basic import (
    IdeogramComponent,
    GencodeComponent,
    TrackComponent,
)

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
    signal_file='path/to/signal.bw',  # Required by TrackComponent
    gencode_annotation_file='path/to/gencode.gtf',  # Required by GencodeComponent
    ideogram_data=read_ideogram('path/to/ideogram.txt'),  # Required by IdeogramComponent (pre-loaded object)
)
```

### Two-Step Workflow (Load then Plot)

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

This is useful when you want to:
- Reuse loaded data for multiple plots with different styling
- Cache data for interactive exploration
- Debug data loading separately from plotting

### Multiple Intervals Per Plot

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
    ideogram_data='path/to/ideogram.txt',
)
```

### Adding Visual Connectors

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

### Parallel Data Loading

For large datasets or slow loaders, enable parallel processing:

```python
data, axes = plotter.plot(
    interval,
    n_cpus=4,  # Use 4 processes for data loading
    **loader_kwargs
)
```

**Note:** Due to multiprocessing overhead, parallel loading may be slower for fast loaders. Best for I/O-bound operations.

## Available Plot Components

### Basic Components

#### TrackComponent
Plots signal tracks from bigWig files.

**Loaders:** `SignalLoader`

**Required arguments:**
- `signal_file` (str): Path to bigWig file

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

#### IdeogramComponent
Displays chromosome ideogram with cytobands.

**Loaders:** `IdeogramLoader`

**Required arguments:**
- `ideogram_data` dict: Result of `genome_tools.plotting.ideogram.read_ideogram` on path to ideogram data file

**Example:**
```python
IdeogramComponent(
    height=0.1,
    margins=(0.1, 0.0),
    ideogram_data=read_ideogram('cytoBand.txt'),
)
```

#### GencodeComponent
Plots gene annotations from GENCODE GTF files.

**Loaders:** `GencodeLoader`

**Required arguments:**
- `gencode_annotation_file` (str): Path to GENCODE GTF file

**Optional plot arguments:**
- `gene_symbol_exclude_regex` (str): Regex to exclude certain gene symbols (default: `r'^ENSG|^MIR|^LINC|.*-AS.*'`)

**Example:**
```python
GencodeComponent(
    height=0.3,
    margins=(0.05, 0.05),
    gencode_annotation_file='gencode.v43.annotation.gtf',
    gene_symbol_exclude_regex=r'^ENSG',  # Exclude Ensembl IDs
)
```

### Sequence Components

#### SequencePlotComponent
Plots DNA sequence as colored nucleotides.

**Loaders:** `FastaLoader`, `OHESequenceLoader`

**Required arguments:**
- `fasta_file` (str): Path to genome FASTA file

**Example:**
```python
from genome_tools.plotting.modular_plot.plot_components.sequence import SequencePlotComponent

SequencePlotComponent(
    height=0.2,
    fasta_file='hg38.fa',
)
```

#### MotifHitsComponent
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

### DHS (DNase Hypersensitivity) Components

#### DHSIndexComponent
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

#### NMFTracksComponent
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

### Advanced: Abstract Base Components

#### SegmentPlotComponent
Abstract base class for plotting genomic segments.

Inherit from this class to create custom segment-based visualizations:

```python
from genome_tools.plotting.modular_plot.plot_components.abstract import SegmentPlotComponent
from genome_tools.plotting.modular_plot import uses_loaders
from genome_tools.plotting.modular_plot.loaders.basic import SegmentsTabixLoader

@uses_loaders(SegmentsTabixLoader)
class CustomSegmentComponent(SegmentPlotComponent):
    __intervals_attr__ = 'intervals'  # Attribute name for segments in DataBundle
    
    def _plot(self, data, ax, **kwargs):
        # Custom plotting logic
        super()._plot(data, ax, color='purple', **kwargs)
        ax.set_ylabel('Custom Segments')
        return ax
```

## Creating Custom Components

### Step 1: Create a Data Loader

Data loaders fetch and process data for plot components. Create a loader by inheriting from `PlotDataLoader`:

```python
from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle

class MyCustomLoader(PlotDataLoader):
    def _load(self, data: DataBundle, my_file: str, threshold: float = 0.5):
        """
        Load and process custom data.
        
        Parameters
        ----------
        data : DataBundle
            The data bundle to populate
        my_file : str
            Path to data file (required)
        threshold : float
            Filtering threshold (optional, default=0.5)
        """
        # Load data based on the interval
        interval = data.interval
        raw_data = load_from_file(my_file, interval)
        
        # Process the data
        processed_data = raw_data[raw_data > threshold]
        
        # Attach to data bundle
        data.custom_data = processed_data
        
        return data
```

**Key points:**
- The `_load` method signature defines required and optional parameters
- Required parameters have no default value
- Optional parameters have default values
- The method must accept `data: DataBundle` as the first argument
- Return the modified `data` object

### Step 2: Create a Simple Loader (No Custom Logic)

For loaders that simply pass through data without processing:

```python
class SimpleLoader(PlotDataLoader):
    required_loader_kwargs = ['my_data', 'my_config']
```

This automatically creates a `_load` method that sets `data.my_data` and `data.my_config`.

### Step 3: Create a Plot Component

Plot components visualize the loaded data:

```python
from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
import matplotlib.pyplot as plt

@uses_loaders(MyCustomLoader)
class MyCustomComponent(IntervalPlotComponent):
    
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax: plt.Axes, color='red', linewidth=2, **kwargs):
        """
        Plot custom data on the provided axes.
        
        Parameters
        ----------
        data : DataBundle
            Data bundle with loaded data
        ax : plt.Axes
            Matplotlib axes to plot on
        color : str
            Line color
        linewidth : float
            Line width
        **kwargs
            Additional matplotlib arguments
        """
        # Access loaded data
        x = data.interval.to_array()
        y = data.custom_data
        
        # Create the plot
        ax.plot(x, y, color=color, linewidth=linewidth, **kwargs)
        ax.set_ylabel('Custom Signal')
        ax.set_ylim(0, max(y) * 1.1)
        
        return ax
```

**Key points:**
- Use `@uses_loaders(...)` decorator to specify required loaders
- Implement `_plot(self, data, ax, **kwargs)` method
- Access loaded data via `data.attribute_name`
- The `@IntervalPlotComponent.set_xlim_interval` decorator automatically sets x-axis limits
- Return the axes object

### Step 4: Use Your Custom Component

```python
from genome_tools.plotting.modular_plot import IntervalPlotter

# Create plotter with custom component
plotter = IntervalPlotter([
    MyCustomComponent(
        height=1.0,
        margins=0.1,
        my_file='data.txt',  # Passed to MyCustomLoader
        threshold=0.7,       # Passed to MyCustomLoader
        color='blue',        # Passed to _plot method
        linewidth=1.5,       # Passed to _plot method
    )
])

# Plot
data, axes = plotter.plot(interval)
```

## Advanced Component Features

### Using Multiple Loaders

Components can use multiple loaders that run sequentially:

```python
@uses_loaders(FastaLoader, OHESequenceLoader, GCContentLoader)
class SequenceWithGCComponent(IntervalPlotComponent):
    def _plot(self, data, ax, **kwargs):
        # Access data from all three loaders
        seq_matrix = data.matrix        # From OHESequenceLoader
        gc_content = data.gc_content    # From GCContentLoader
        
        # Plot sequence with GC content overlay
        seq_plot(seq_matrix, ax=ax)
        ax_gc = ax.twinx()
        ax_gc.plot(gc_content, color='green', alpha=0.5)
        
        return ax
```

### Overriding Loader Parameters

Parameters can be overridden at multiple levels (from lowest to highest priority):

1. **Loader defaults** - Defined in the loader's `_load` method signature
2. **Global overrides** - Passed to `IntervalPlotter.__init__`
3. **Component-level overrides** - Passed to component `__init__`
4. **Loader-specific component overrides** - Using loader name as key
5. **Runtime overrides** - Passed to `plot()` or `get_interval_data()`
6. **Runtime loader-specific overrides** - Using loader name as key

```python
# Global override (applies to all components)
plotter = IntervalPlotter(
    [comp1, comp2],
    signal_file='default.bw',  # Used by both components
)

# Component-level override
comp = TrackComponent(
    signal_file='specific.bw',  # Overrides global default
)

# Loader-specific override (when loader args conflict)
from genome_tools.plotting.modular_plot.loaders.basic import SignalLoader

comp = MyComponent(
    threshold=0.5,  # Used by multiple loaders
    SignalLoader={'threshold': 0.8},  # Only for SignalLoader
)

# Runtime override
data, axes = plotter.plot(
    interval,
    signal_file='runtime.bw',  # Highest priority
)
```

### Creating Component Variants with `with_loaders`

Create new component classes with different loader combinations:

```python
from genome_tools.plotting.modular_plot.plot_components.basic import TrackComponent
from genome_tools.plotting.modular_plot.loaders.basic import AverageSignalLoader

# Create a variant that averages multiple signal files
AverageTrackComponent = TrackComponent.with_loaders(
    AverageSignalLoader,
    new_class_name='AverageTrackComponent'
)

# Use the new component
comp = AverageTrackComponent(
    height=1.5,
    signal_files=['file1.bw', 'file2.bw', 'file3.bw'],  # Required by AverageSignalLoader
    smooth=True,
    bandwidth=150,
)
```

### Custom Component with Configurable Loader

```python
class FlexibleComponent(IntervalPlotComponent):
    """A component that can use different loaders based on configuration."""
    pass

# Create variants for different data sources
BigWigComponent = FlexibleComponent.with_loaders(
    SignalLoader,
    new_class_name='BigWigComponent'
)

BedGraphComponent = FlexibleComponent.with_loaders(
    BedGraphLoader,
    new_class_name='BedGraphComponent'
)
```

### Accessing Component Axes

The `plot()` and `plot_interval()` methods return a named tuple of axes:

```python
data, axes = plotter.plot(interval)

# Access by component name
ideogram_ax = axes.IdeogramComponent
gencode_ax = axes.GencodeComponent
track_ax = axes.TrackComponent

# Access by index
first_ax = axes[0]
second_ax = axes[1]

# Iterate
for ax in axes:
    ax.set_xlabel('Genomic Position')
```

### Component Naming

By default, components are named after their class. Customize names to use multiple instances:

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

### Debugging Data Loading

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

### Logging

Enable logging to debug component and loader behavior:

```python
import logging

# Component with debug logging
comp = TrackComponent(
    logger_level=logging.DEBUG,
    height=1.5,
)

# Plotter with info logging
plotter = IntervalPlotter(
    [comp1, comp2],
    logger_level=logging.INFO,
)
```

## Best Practices

### Loader Design
- Keep loaders focused on a single data source or transformation
- Make parameters explicit in the `_load` signature
- Add type hints and docstrings
- Validate input data and provide helpful error messages
- Avoid side effects; always return the modified `data` object

### Component Design
- Separate data loading (loaders) from visualization (components)
- Use `@IntervalPlotComponent.set_xlim_interval` for x-axis consistency
- Accept matplotlib styling arguments as `**kwargs`
- Return the axes object from `_plot`
- Set reasonable default heights and margins

### Performance
- Use `n_cpus > 1` only for I/O-bound loaders (file reading)
- Cache loaded data for multiple plots with different styles
- Use the two-step workflow (`get_interval_data` + `plot_interval`) for iteration

### Organization
- Group related loaders in a module (e.g., `loaders/chipseq.py`)
- Group related components in a module (e.g., `plot_components/epigenomics.py`)
- Use abstract base components for common patterns
- Document required loader arguments in component docstrings

## Complete Example: Custom Epigenomics Component

```python
from genome_tools.plotting.modular_plot import PlotDataLoader, IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.utils import DataBundle
from genome_tools.data.extractors import BigwigExtractor
import numpy as np
import matplotlib.pyplot as plt

# Loader: Fetch ChIP-seq data
class ChIPSeqLoader(PlotDataLoader):
    def _load(self, data: DataBundle, chip_file: str, input_file: str = None, log_scale: bool = False):
        """
        Load ChIP-seq signal and optionally subtract input.
        """
        with BigwigExtractor(chip_file) as extractor:
            chip_signal = np.nan_to_num(extractor[data.interval])
        
        if input_file:
            with BigwigExtractor(input_file) as extractor:
                input_signal = np.nan_to_num(extractor[data.interval])
            signal = chip_signal - input_signal
        else:
            signal = chip_signal
        
        if log_scale:
            signal = np.log1p(signal)
        
        data.chip_signal = signal
        return data

# Component: Visualize ChIP-seq data
@uses_loaders(ChIPSeqLoader)
class ChIPSeqComponent(IntervalPlotComponent):
    
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, color='darkblue', fill=True, **kwargs):
        """
        Plot ChIP-seq signal as a line or filled area.
        """
        x = data.interval.to_array()
        y = data.chip_signal
        
        if fill:
            ax.fill_between(x, y, color=color, alpha=0.5, **kwargs)
        else:
            ax.plot(x, y, color=color, **kwargs)
        
        ax.set_ylabel('ChIP-seq Signal')
        ax.set_ylim(0, np.percentile(y, 99) * 1.1)
        ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
        
        return ax

# Usage
from genome_tools import GenomicInterval
from genome_tools.plotting.modular_plot import IntervalPlotter

interval = GenomicInterval('chr1', 1_000_000, 1_100_000)

plotter = IntervalPlotter([
    ChIPSeqComponent(
        height=1.5,
        margins=(0.1, 0.1),
        chip_file='h3k27ac.bw',
        input_file='input.bw',
        log_scale=True,
        color='red',
        fill=True,
    )
])

data, axes = plotter.plot(interval)
plt.savefig('chipseq_track.pdf')
```

## Summary

The modular plot framework enables:
- **Composability** - Mix and match components without code changes
- **Reusability** - Share loaders across different plot types
- **Flexibility** - Override parameters at multiple levels
- **Scalability** - Parallel data loading for complex plots
- **Maintainability** - Clean separation of data loading and visualization

Start with the built-in components, then create custom components as your visualization needs grow. The framework handles the complexity of figure layout, data management, and parameter passing, letting you focus on the science

