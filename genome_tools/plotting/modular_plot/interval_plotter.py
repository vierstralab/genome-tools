from typing import Sequence, List
from collections import namedtuple
from matplotlib import gridspec
from matplotlib import pyplot as plt
import inspect

from genome_tools import GenomicInterval
import pandas as pd

from genome_tools.plotting.modular_plot.api import PlotComponent, IntervalPlotComponent
from genome_tools.plotting.modular_plot.utils import LoggerMixin, DataBundle
from genome_tools.plotting.connectors import connect_axes_lines, connect_axes_area
from genome_tools.plotting.utils import format_axes_to_interval


class PlotComponentManager(LoggerMixin):
    """
    Base class to manage plot components. Ensures unique component names and stores their order
    """
    def __init__(self, plot_components: Sequence[PlotComponent], logger_level=None):
        super().__init__(logger_level=logger_level)
        self.component_names = [c.name for c in plot_components]

        repeating_names = self._check_unique_component_names(self.component_names)
        if len(repeating_names) > 0:
            message = f"""Component names must be unique. 
            Please, explicitly set the name for each repeating component class. 
            Repeating component names: {repeating_names}"""
            self.logger.error(message)
            raise ValueError(message)

        self.CompTuple = namedtuple('ComponentNamesTuple', self.component_names)

        self.plot_components = self.CompTuple(*plot_components)

        self._validate_overlapping_component_kwargs()
    
    def _validate_overlapping_component_kwargs(self):
        """
        Validate that there are no overlapping kwargs between components.
        """
        kwarg_to_components = {}

        for component in self.plot_components:
            component: PlotComponent
            comp_name = component.name
            kwarg_set = component.__loader_kwargs_signature__

            for kw in kwarg_set:
                kwarg_to_components.setdefault(kw, []).append(comp_name)
        
        overlapping = {
            kw: comps for kw, comps in kwarg_to_components.items()
            if len(comps) > 1
        }

        if len(overlapping) > 0:
            msg_lines = ["Found overlapping component loader kwargs:"]
            for kw, comps in overlapping.items():
                msg_lines.append(f"Argument '{kw}' used by: {', '.join(comps)}")
            
            full_msg = "\n".join(msg_lines)
            self.logger.warning(full_msg)

    @staticmethod
    def _check_unique_component_names(names):
        """
        Check if the component names are unique. And return the non-unique names.
        """
        value_counts = pd.Series(names).value_counts()
        return value_counts[value_counts > 1].index.tolist()

    @staticmethod
    def _convert_component_to_name(component):
        """
        Convert a vertical plot component to its name.
        """
        if isinstance(component, PlotComponent):
            return component.name
        elif isinstance(component, type) and issubclass(component, PlotComponent):
            return component.__name__
        else:
            assert isinstance(component, str)
            return component
        
    def _sort_components(self, components: Sequence[str]):
        """
        Sort the components in the order of self.plot_components.
        """
        if not all(c in self.component_names for c in components):
            self.logger.error(f"Invalid component names. Must be one of the following: {self.component_names}")
            raise AssertionError
        return [c for c in self.component_names if c in components]
    
    def _get_interim_components(self, component1, component2):
        idx1 = self.component_names.index(component1)
        idx2 = self.component_names.index(component2)
        return self.component_names[idx1:idx2 + 1]

    
    def get_component_axes(self, component_axes, component):
        """
        Get the axes object for a specific component.
        If the component's plot method returns multiple axes, return the last one.
        """
        try:
            axes = getattr(component_axes, component)
            if not isinstance(axes, plt.Axes):
                # workaround when PlotComponent returns multiple axes
                assert len(axes) > 0
                assert isinstance(axes[0], plt.Axes)
                self.logger.warning(f"PlotComponent {component} returned multiple axes. Using the first one when adding a connector.")
                return axes[0]
            return axes
        except AttributeError:
            self.logger.error(f"Component axes do not have attribute {component}.")
            raise

class VerticalConnectorMixin(PlotComponentManager):
    def get_connecting_components(self, start_component, end_component):
        if start_component is None:
            start_component = self.component_names[0]
        if end_component is None:
            end_component = self.component_names[-1]
        
        if start_component == end_component:
            self.logger.error("Start and end components cannot be the same.")
            raise ValueError("Start and end components cannot be the same.")

        components = [
            self._convert_component_to_name(c) for c in 
            [start_component, end_component]
        ]
        components = self._sort_components(components)
        components = self._get_interim_components(components[0], components[1])
        return components

    def plot_connector(
            self,
            component_axes,
            positions: List[int]=None,
            interval: GenomicInterval=None,
            start_component=None,
            end_component=None,
            type='area',
            extend_to_top=False,
            extend_to_bottom=False,
            **kwargs,
        ):
        """
        Add a connector between two vertical plot components. Can run multiple times
        
        Parameters
        ----------
        start_component : Union[str, VerticalPlotComponent], optional
            The starting vertical plot component for the connector. If None uses the first component.
        
        end_component : Union[str, VerticalPlotComponent], optional
            The ending vertical plot component for the connector. If None uses the last component.
        
        type : str, optional
            The type of connector. Must be one of 'line' or 'area'.
            Default is 'area'.
            'line' is used to connect single points across the components with lines.
            'area' is used to connect intervals across the components with shaded areas.

        extend_to_top : bool, optional
            If True, the connector extends to the top of the first component axes.
            Default is False.

        extend_to_bottom : bool, optional
            If True, the connector extends to the bottom of the last component axes.
            Default is False.

        **kwargs
            Additional keyword arguments to pass to the connector plot method.
            See VerticalAxesConnector for more details.

        Returns
        -------
        connectors : list[VerticalAxesConnector]
            The connectors that were added.

        Usage
        -----
        # Add connectors between components
        connectors = connector_plotter.add_connectors('ComponentA', 'ComponentB', type='line', x=[1, 2, 3])
        """
        
        components = self.get_connecting_components(start_component, end_component)
        n = len(components)

        location_kwargs = self.parse_kwargs(positions, interval, type)
        for i in range(n - 1):
            component1 = components[i]
            component2 = components[i + 1]
            connector_kwargs = {
                **kwargs,
                **location_kwargs,
                'extend_to_top': extend_to_top if i == 0 else False,
                'extend_to_bottom': extend_to_bottom if i == n - 2 else True,
            }
            ax1 = self.get_component_axes(component_axes, component1)
            ax2 = self.get_component_axes(component_axes, component2)
            # TODO: resolve component axes here instead of VerticalAxesConnector
            if type == 'area':
                connect_axes_area(ax_top=ax1, ax_bottom=ax2, **connector_kwargs)
            elif type == 'line':
                connect_axes_lines(ax_top=ax1, ax_bottom=ax2, **connector_kwargs)

    
    @staticmethod
    def parse_kwargs(positions, interval, type):
        assert type in ['line', 'area'], "Connector type must be one of 'line' or 'area'."
        kwargs = {}
        if positions is None:
            if interval is not None:
                if not isinstance(interval, GenomicInterval):
                    # add key parsing for the interval
                    raise ValueError("Interval must be a GenomicInterval object.")
                positions = [interval.start, interval.end]
            else:
                raise ValueError("Either positions or interval must be provided.") 
        if type == "area":
            assert len(positions) == 2, "For area connectors, positions must be a list of two elements: [start, end]."
            kwargs['x1'] = positions[0]
            kwargs['x2'] = positions[1]
        elif type == "line":
            kwargs['x'] = positions
        return kwargs

class IntervalPlotter(VerticalConnectorMixin):
    """
    Class to plot a genomic interval with multiple vertical plot components.

    Parameters
    ----------
    plot_components : Sequence[VerticalPlotComponent]
        The vertical plot components to plot in the interval.

    inches_per_unit : float, optional
        The number of inches per unit for the figure size.
        Default is 1.0.

    width : float, optional
        The width of the figure in inches * inches_per_unit.
        Default is 2.5.
    
    **data_kwargs : dict
        Keyword arguments to pass to the loaders.

    Usage
    -----

    # Define plot components
    plot_components = [
        IdeogramComponent(height=0.1, margins=(0.1, 0.1), interval_key='gene'),
        GencodeComponent(height=0.2, margins=(0.1, 0.1), interval_key='gene'),
        FinemapComponent(height=0.2, margins=(0.1, 0.1), interval_key='dhs'),
        DNaseTracksComponent(height=1.5, margins=(0.1, 0.1), interval_key='dhs', component_data=component_data),
        DHSIndexComponent(height=0.1, margins=(0.1, 0.0), interval_key='dhs'),
        DHSLoadingsComponent(height=0.2, margins=(0.0, 0.1), interval_key='dhs', component_data=component_data),
        FootprintsComponent(height=0.1, margins=(0.1, 0.1), interval_key='footprint'),
        MotifComponent(height=0.2, margins=(0.1, 0.1), interval_key='footprint'),
        CAVComponent(height=0.2, margins=(0.1, 0.1), interval_key='footprint'),
    ]

    # Create an interval plotter
    interval_plotter = IntervalPlotter(plot_components)

    # plot with just one command
    component_axes, data, connectors = interval_plotter.plot(interval)

    # Alternatively you can run 2 separate commands for more control:
    # Get the data for the interval
    data = interval_plotter.get_interval_data(interval)

    # Plot the interval
    component_axes = interval_plotter.plot_interval(data)
    """
    def __init__(self, plot_components: Sequence[IntervalPlotComponent],
                inches_per_unit=1.0, width=2.5, **data_kwargs):
        super().__init__(plot_components)
        self.data_kwargs = data_kwargs

        self.inches_per_unit = inches_per_unit
        self.width = width
        self.gridspec = self.setup_default_gridspec()
    
    def __repr__(self):
        return f"IntervalPlotter(Plotcomponents=[{', '.join(self.component_names)}])"
    
    def setup_default_gridspec(self):
        """
        Setup the GridSpec for the vertical components in the figure.
        Components are plotted in a vertical stack with specified heights and margins.
        """
        height_ratios = [
            x
            for c in self.plot_components
            for x in [c.margin_top, c.height, c.margin_bottom]
        ]
        return gridspec.GridSpec(len(height_ratios), 1, height_ratios=height_ratios, hspace=0)


    def get_all_component_gridspecs(self):
        """
        Get the GridSpecs for each vertical plot component.
        """
        return self.CompTuple(*[self.gridspec[3 * i + 1, :] for i in range(len(self.plot_components))])
    
    def get_gridspec_for_component(self, component, include_top_margin=False,
                                   include_bottom_margin=False):
        """
        Get the GridSpec slice for a specific vertical plot component.
        """
        index = self.component_names.index(self._convert_component_to_name(component))
        start = 3 * index if include_top_margin else 3 * index + 1
        end = 3 * index + 3 if include_bottom_margin else 3 * index + 2
        return self.gridspec[start:end, :]
    
    def setup_default_figure(self):
        """
        Setup a default figure with the appropriate size for the vertical components.
        """
        return plt.figure(
            figsize=(self.inches_per_unit * self.width,
                sum(
                    x for c in self.plot_components 
                    for x in [c.margin_top, c.height, c.margin_bottom]
                ) * self.inches_per_unit
            )
        )

    def get_interval_data(
            self,
            interval: GenomicInterval,
            **data_kwargs
        ):
        """
        Get the data for the specified interval(s) and plot components.

        Parameters
        ----------
        interval : GenomicInterval or dict {key: GenomicInterval}
            The genomic interval to plot.
            If a dict is provided, the component-specific interval key is used to extract the interval.

        plot_components : Sequence[VerticalPlotComponent]
            The vertical plot components to plot.

        **data_kwargs : dict
            Keyword arguments to pass to the loaders function.

        Returns
        -------
        data : Iter[DataBundle]
            A Iter of DataBundle objects containing the data for each plot component
        """
        common_kwargs = set(data_kwargs) & set(self.data_kwargs)
        if common_kwargs:
            self.logger.debug(
                f"Found {len(common_kwargs)} overlapping data kwargs: {list(common_kwargs)}"
            )
            self.logger.debug("Using values passed to get_interval_data function.")

        result = []
        # Use multiprocessing or threading here if needed in the future
        for component in self.plot_components:
            component: PlotComponent
            component_interval = self._parse_interval(
                interval,
                getattr(component, 'interval_key', None)
            )
            if component_interval is None:
                self.logger.error(f"Component {component.name} requires an interval_key to extract the interval from the provided dict.")
                raise ValueError(f"Component {component.name} requires an interval_key to extract the interval from the provided dict.")

            data = DataBundle(
                interval=component_interval
            )
            result.append(
                component.load_data(
                    data,
                    **{**self.data_kwargs, **data_kwargs}
                )
            )
        return self.CompTuple(*result)

    def plot_interval(self, data: Sequence[DataBundle], fig=None, gridspecs=None, **kwargs):
        """
        Plot the genomic interval with all the provided vertical plot components.

        Parameters
        ----------
        data : Sequence[DataBundle] # Named tuple
            The data bundles for each vertical plot component.
            Should be a named tuple with the same names as the component names.
            Generated by get_interval_data.
        
        fig : Figure, optional
            The figure to plot the interval.
            If None, a new figure is created.

        gridspecs : Sequence[GridSpec], optional
            The GridSpecs for each vertical plot component.
            If None, the default GridSpecs are used.

        Usage
        -----
        # Get the data for the interval
        data = interval_plotter.get_interval_data(interval, plot_components)
        # Plot the interval
        component_axes = interval_plotter.plot_interval(data)
        """
        if fig is None:
            fig = self.setup_default_figure()
        
        component_axes = []

        if gridspecs is None:
            gridspecs = self.get_all_component_gridspecs()

        assert len(gridspecs) == len(self.plot_components)

        try:
            filtered_data = [getattr(data, c) for c in self.component_names]
        except AttributeError:
            missing_components = [c for c in self.component_names if not hasattr(data, c)]
            self.logger.error("Provided data does not have the correct attributes for the following components: " + ", ".join(missing_components))
            raise

        for gs, component, data_bundle in zip(gridspecs, self.plot_components, filtered_data):
            component: IntervalPlotComponent
            ax = fig.add_subplot(gs)
            format_axes_to_interval(ax, data_bundle.interval, axis='x')
            component_axes.append(component.plot(data_bundle, ax=ax, **kwargs))
        component_axes = self.CompTuple(*component_axes)

        return component_axes
    
    def plot(self, interval: GenomicInterval, fig=None, gridspecs=None, data_kwargs=None, **plot_kwargs):
        """
        Plot the genomic interval with all the provided vertical plot components.

        Parameters
        ----------
        interval : GenomicInterval or dict {key: GenomicInterval}
            The genomic interval to plot.
            If a dict is provided, the component-specific interval key is used to extract the interval.

        fig : Figure, optional
            The figure to plot the interval.
            If None, a new figure is created.

        gridspecs : Sequence[GridSpec], optional
            The GridSpecs for each vertical plot component.
            If None, the default GridSpecs are used.

        **data_kwargs : dict
            Keyword arguments to pass to the loaders function.

        Returns
        -------
        component_axes : Iter[Axes]
            The axes for each vertical plot component.

        data : Iter[DataBundle]
            The data bundles for each vertical plot component.
        
        Usage
        -----
        # Plot the interval
        component_axes, data = interval_plotter.plot(interval)
        """
        if data_kwargs is None:
            data_kwargs = {}
        
        data = self.get_interval_data(interval, **data_kwargs)
        component_axes = self.plot_interval(data, fig=fig, gridspecs=gridspecs, **plot_kwargs)
        return data, component_axes

    @staticmethod
    def _parse_interval(interval, interval_key: str):
        """
        Parse the interval argument to a GenomicInterval object.
        If a dict is provided, the interval_key is used to extract the interval.
        """
        if isinstance(interval, dict):
            try:
                return interval[interval_key]
            except KeyError:
                raise ValueError(f"Interval key '{interval_key}' not found.")
        else:
            return interval