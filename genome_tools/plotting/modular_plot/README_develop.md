# Develop
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

# Plotter with info logging
plotter = IntervalPlotter(
    [comp1, comp2],
    logger_level=logging.INFO,
)
```

Per component returned axes match with whatever is returned by `_plot`. When returning multiple axes, first axis should be the input ax to ensure compatibility with connectors.