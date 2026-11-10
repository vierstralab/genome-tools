# Develop

## Table of contents
- [Creating Custom Components](#creating-custom-components)
  - [Step 1: Create a Data Loader](#step-1-create-a-data-loader)
  - [Step 2: Create a Simple Loader (No Custom Logic)](#step-2-create-a-simple-loader-no-custom-logic)
  - [Step 3: Create a Plot Component](#step-3-create-a-plot-component)
  - [Step 4: Use Your Custom Component](#step-4-use-your-custom-component)
  - [Example: heatmap of several BigWigs](#example-heatmap-of-several-bigwigs)
- [Advanced: Abstract Base Components](#advanced-abstract-base-components)
  - [SegmentPlotComponent](#segmentplotcomponent)
- [Advanced Component Features](#advanced-component-features)
  - [Using Multiple Loaders](#using-multiple-loaders)
  - [Creating Component Variants with `with_loaders`](#creating-component-variants-with-with_loaders)
  - [Custom Component with Configurable Loader](#custom-component-with-configurable-loader)
- [Best Practices](#best-practices)
  - [Loader Design](#loader-design)
  - [Component Design](#component-design)
  - [Logging](#logging)

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
import numpy as np
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
        x = np.arange(data.interval.start, data.interval.end)
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

### Example: heatmap of several BigWigs

A complete component that runs on public data: DNase-seq density of three samples as a heatmap, one row per sample, under the gene annotation.

```python
import numpy as np
from genome_tools import GenomicInterval
from genome_tools.data.extractors import BigwigExtractor
from genome_tools.plotting import cm2inch
from genome_tools.plotting.modular_plot import IntervalPlotter, IntervalPlotComponent, PlotDataLoader, uses_loaders
from genome_tools.plotting.modular_plot.plot_components.basic import GencodeComponent


class MultiSignalLoader(PlotDataLoader):
    """Read per-base signal from several BigWigs (`signal_files`), averaged in `bin_size` bp bins.
    Sets `data.matrix` (n_files, n_bins) and `data.row_labels`."""
    def _load(self, data, signal_files, row_labels=None, bin_size=1):
        rows = []
        for f in signal_files:
            with BigwigExtractor(f) as bw:
                rows.append(np.nan_to_num(bw[data.interval]))
        matrix = np.stack(rows)
        if bin_size > 1:
            n_bins = matrix.shape[1] // bin_size
            matrix = matrix[:, :n_bins * bin_size].reshape(len(rows), n_bins, bin_size).mean(axis=2)
        data.matrix = matrix
        data.row_labels = row_labels if row_labels is not None else [str(i) for i in range(len(rows))]
        return data


@uses_loaders(MultiSignalLoader)
class HeatmapComponent(IntervalPlotComponent):
    """Several signals over an interval as a heatmap, one row per file."""
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, cmap='Reds', **kwargs):
        ax.imshow(
            data.matrix,
            aspect='auto',
            interpolation='nearest',
            cmap=cmap,
            extent=(data.interval.start, data.interval.end, len(data.matrix), 0),  # genomic x coordinates
            **kwargs,
        )
        ax.set_yticks(np.arange(len(data.row_labels)) + 0.5)
        ax.set_yticklabels(data.row_labels)
        return ax


samples = ['AG95384', 'AG95451', 'AG80853']
plotter = IntervalPlotter(
    [
        GencodeComponent(height=cm2inch(0.6)),
        HeatmapComponent(
            signal_files=[f'https://regulome.altius.org/data/per_sample/{s}/{s}.normalized_density.bw' for s in samples],
            row_labels=samples, bin_size=50,                # loader args
            cmap='magma_r', vmin=0, vmax=10,                # plot kwargs, passed to ax.imshow
            height=cm2inch(1),
        ),
    ],
    gencode_annotation_file='https://resources.altius.org/~sabramov/collaborations/vinson/examples/gencode.v49.basic.annotation.sorted.gff3.gz',
    width=cm2inch(8),
)
data, axes = plotter.plot(GenomicInterval.from_ucsc('chr2:60,485,000-60,505,000'))
```


## Advanced: Abstract Base Components
You can inherit from these components to make basic plots

### SegmentPlotComponent
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

### Creating Component Variants with `with_loaders`

Create new component classes with different loader combinations:

```python
from genome_tools.plotting.modular_plot.plot_components.basic import TrackComponent
from genome_tools.plotting.modular_plot.loaders.basic import AverageSignalLoader

# Create a variant that averages multiple signal files
AverageTrackComponent = TrackComponent.with_loaders(
    AverageSignalLoader,
    new_class_name='AverageTrackComponent',
    doc="Average of several signal files (`signal_files`), optionally smoothed.",  # docstring of the new class
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

### Logging

Enable logging to debug component and loader behavior:

```python
import logging

# Component with debug logging
comp = TrackComponent(
    logger_level=logging.DEBUG,
    height=1.5,
)

# Plotter with info logging: extra IntervalPlotter kwargs go to the loaders,
# so set the level on its logger instead
plotter = IntervalPlotter([comp])
plotter.logger.setLevel(logging.INFO)
```

Per component returned axes match with whatever is returned by `_plot`. When returning multiple axes, first axis should be the input ax to ensure compatibility with connectors.