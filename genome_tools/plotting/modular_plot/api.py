import inspect

from matplotlib import pyplot as plt
from typing import List, Type, Union, Sequence
import types
from makefun import wraps, add_signature_parameters, remove_signature_parameters, with_signature
import warnings


from genome_tools.plotting.utils import clear_spines, format_axes_to_interval
from genome_tools import GenomicInterval

from genome_tools.plotting.modular_plot.utils import LoggerMixin, DataBundle, RequiredArgument


class PlotDataLoader(LoggerMixin):
    """
    Base class for all data loaders.
    Each loader should implement the `_load` method to load and filter data based on the interval.

    By default (if not implemented), the _load method is created from the required_loader_kwargs. The generated _load method takes required_loader_kwargs as arguments and sets them as fields of the data object.
    """
    required_loader_kwargs = []

    def __init__(self, logger_level=None):
        LoggerMixin.__init__(self, logger_level=logger_level)

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls._load is PlotDataLoader._load:
            cls._set_default_load()
        elif len(cls.required_loader_kwargs) > 0:
            warnings.warn(
                f"Both required_loader_kwargs and _load are specified for loader {cls.__name__}. required_loader_kwargs are ignored."
            )

    @classmethod
    def _set_default_load(cls):
        """
        Set the default load method if not already set.
        """
        if len(cls.required_loader_kwargs) == 0:
            warnings.warn(
                f"Loader {cls.__name__} has no required_loader_kwargs and no _load method implemented. The loader does not modify the data."
            )

        args_string = ', '.join(cls.required_loader_kwargs)
        @with_signature(f"_load(self, data, {args_string})")
        def _default_load(self, data: DataBundle, **kwargs):
            """
            Default load method when _load is not implemented.

            Sets the fields of the data object based on the required_loader_kwargs.
            """
            for field, value in kwargs.items():
                setattr(data, field, value)
            return data

        cls._load = types.MethodType(_default_load, cls)

    @classmethod
    def get_fullargspec(cls):
        """
        Get the fullargspec of the load method.
        Returns a dictionary of the arguments and their default values,
        excluding the 'self' and 'data' arguments.
        """
        fullargspec = inspect.getfullargspec(cls._load)
        if fullargspec.varkw is not None or fullargspec.varargs is not None:
            raise ValueError(f"{cls.__name__} '_load' method should not have *args or **kwargs.")
        loader_args = fullargspec.args
        loader_defaults = fullargspec.defaults or []
        extended_defaults = [RequiredArgument()] * (len(loader_args) - len(loader_defaults)) + list(loader_defaults)
        signature = dict(zip(loader_args, extended_defaults))
        if 'data' not in signature:
            raise ValueError(f"{cls.__name__} '_load' method should have 'data' argument.")
        return {k: v for k, v in signature.items() if k not in ['self', 'data']}

    def _load(self, data: DataBundle):
        """
        This class should load the data
        """
        raise NotImplementedError
    
    def load(self, data: DataBundle, **loader_kwargs):
        """
        Validate the preprocessor fields and call the load method.
        """
        self._validate(**loader_kwargs)
        loaded_data = self._load(data, **loader_kwargs)
        if loaded_data is None:
            loaded_data = data
        return loaded_data

    def _validate(self, **loader_kwargs):
        """
        Validate that the preprocessor has all the required fields.
        """
        missing_args = [
            arg for arg in loader_kwargs 
            if isinstance(loader_kwargs[arg], RequiredArgument)
        ]
        if missing_args:
            raise ValueError(f"Loader {self.__class__.__name__} is missing required argument(s): {', '.join(missing_args)}")


class PlotComponent(LoggerMixin):
    """
    Base class for all plot components.
    Each plot component should implement the `_plot` method,
    which takes the data and an axis object and plots the data on it.
    Use the `uses_loaders` decorator to specify which data loaders are required by the plot component or use the `with_loaders` class method to create a new class with the specified loaders.
    """
    __loader_kwargs_signature__ = {} # contains all loader kwargs for updated signature
    __required_loaders__: List[Type[PlotDataLoader]] = []
    __original_init__ = None

    def __init__(self, name=None, logger_level=None, **kwargs):
        LoggerMixin.__init__(self, logger_level=logger_level)

        if name is None:
            name = self.__class__.__name__
        self.name = name

        # Separate loader kwargs from plot kwargs
        self.loader_kwargs = {
            key: kwargs.pop(key, v)
            for key, v in self.__class__.__loader_kwargs_signature__.items()
        }

        self.plot_kwargs = kwargs

    @classmethod
    def with_loaders(cls, *loaders, new_class_name=None):
        """
        Create a new class that inherits from the current class
        but requires the specified loaders.
        """
        if new_class_name is None:
            new_class_name = cls.__name__
        new_class = type(
            new_class_name,
            (cls,),
            {}
        )
        return uses_loaders(*loaders)(new_class)
    
    def load_data(self, data: DataBundle, **loader_kwargs):
        """
        Modifies data for the plot component using the required loaders load function.
        """
        for LoaderClass in self.__required_loaders__:
            all_loader_kwargs = {**self.loader_kwargs, **loader_kwargs}
            all_loader_kwargs = {
                k: v for k, v in all_loader_kwargs.items()
                if k in LoaderClass.get_fullargspec()
            }

            loader = LoaderClass(
                logger_level=self.logger.level
            )
            data = loader.load(data, **all_loader_kwargs)
            data.processed_loaders.append(LoaderClass)
        return data

    def plot(self, data, ax, **plot_kwargs):
        """
        Wrapper for the plot method to pass the plot_kws to the plot method.
        Also supports axes set methods e.g. xlim -> ax.set_xlim
        """
        kws = {**self.plot_kwargs, **plot_kwargs}

        set_methods = [method[len('set_'):] for method in dir(plt.Axes) if method.startswith('set_')]
        axes_kws = kws.pop('axes_kwargs', {})
        axes_kws = {key: axes_kws[key] for key in axes_kws if key in set_methods}

        axes = self._plot(data, ax, **kws)
        if axes is None:
            axes = ax
        for key in axes_kws:
            getattr(ax, 'set_' + key)(axes_kws[key])

        return axes

    def _plot(self, data, ax, **kwargs):
        """
        Abstract plot method to be implemented by specific plot components.

        Should not include any axes set methods as arguments,
        as they will be intercepted by the plot method.
        E.g. xlim, ylim, xlabel, etc.
        """
        raise NotImplementedError("Plot method should be implemented in subclasses.")


class IntervalPlotComponent(PlotComponent):
    """
    An extension of the PlotComponent class that plots data vertically.

    Accepts height and margins arguments to control the size of the plot,
    that are only used when plotting with the IntervalPlotter.plot_interval method.
    Additionaly, the interval_key argument can be used to specify different
    intervals for each plot component.
    Name field is used by the IntervalPlotter to refer to the plot component,
    defaulting to the class name.
    """
    def __init__(self, 
                 height: float = 1.0, 
                 margins: Union[float, Sequence[float]] = 0.1, 
                 interval_key=None,
                 **kwargs):
        self.height = height
        self.margin_top, self.margin_bottom = self._parse_margins(margins)
        self.interval_key = interval_key

        super().__init__(**kwargs)


    def _parse_margins(self, margins):
        """Helper function to parse margins input."""
        if isinstance(margins, Sequence) and not isinstance(margins, str):
            if len(margins) == 2:
                return margins
            else:
                raise ValueError("If margins is a sequence, it must have exactly two elements.")
        return margins, margins


    @staticmethod
    def set_xlim_interval(func):
        """
        Decorator to set the x-axis limits for the plot.
        The interval attribute of the data object is used to set the limits.
        """
        def wrapper(self, data, ax, **kwargs):
            # TODO: figure out ho to store data with named arguments
            # if not hasattr(data, 'interval') or not data.interval:
            #     raise ValueError("Data must have 'interval' attribute.")
            format_axes_to_interval(ax, data.interval)
            return func(self, data, ax, **kwargs)
        return wrapper
    
    def add_axes_at_summits(
        self, 
        genomic_intervals: List[GenomicInterval],
        interval: GenomicInterval,
        bp_width=150,
        summit_field='dhs_summit',
        ax=None
    ):
        """
        Add axes at the summit points of the genomic intervals.

        Parameters
        ----------
        genomic_intervals : list of GenomicInterval
            List of genomic intervals.
        interval : GenomicInterval
            Interval to restrict the axes.
        bp_width : float, optional
            Width in base pairs around the summit to plot.
        summit_field : str, optional
            Field of the genominc_intervals to use as the summit to plot around. Default is 'dhs_summit'.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, use the current axes.

        Returns
        -------
        axes : list of matplotlib.axes.Axes
            List of axes added at the summit points of the genomic intervals.
        """
        modified_genomic_intervals = []
        for genomic_interval in genomic_intervals:
            try:
                summit = getattr(genomic_interval, summit_field)
            except AttributeError:
                summit = (genomic_interval.start + genomic_interval.end) / 2

            modified_genomic_intervals.append(
                GenomicInterval(
                    genomic_interval.chrom,
                    summit,
                    summit + 1
                ).widen((bp_width - 1) // 2)
            )
        return self.add_axes_at_intervals(
            modified_genomic_intervals,
            interval,
            ax=ax
        )

    @staticmethod
    def add_axes_at_intervals(
        genomic_intervals: List[GenomicInterval],
        interval: GenomicInterval,
        ax=None
    ):
        """
        Add axes at the middle points of the genomic intervals.

        Parameters
        ----------
        genomic_intervals : list of GenomicInterval
            List of genomic intervals.
        interval : GenomicInterval
            Interval to restrict the axes.
        bp_width : float, optional
            Width in base pairs around the summit to plot. If None, plot the whole genomic interval.
        summit_field : str, optional
            Field of the genominc_intervals to use as the summit
            to plot around. Default is 'dhs_summit'. Only used if bp_width is not None.
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, use the current axes.

        Returns
        -------
        axes : list of matplotlib.axes.Axes
            List of axes added at the middle points of the genomic intervals.
        """
        if ax is None:
            ax = plt.gca()

        axes = []
        for genomic_interval in genomic_intervals:
            parent_pos = ax.get_position()
            # work in start coordinates
            x0_rel = (genomic_interval.start - interval.start) / (interval.end - 1 - interval.start)
            x1_rel = (genomic_interval.end - 1 - interval.start) / (interval.end - 1- interval.start)

            new_axes_width = parent_pos.width * (x1_rel - x0_rel)
            new_axes_height = parent_pos.height
            new_axes_x = parent_pos.x0 + (x0_rel * parent_pos.width)
            new_axes_y = parent_pos.y0

            new_ax = ax.get_figure().add_axes([new_axes_x, new_axes_y, new_axes_width, new_axes_height])
            new_ax.set_xticks([])
            new_ax.set_yticks([])
            clear_spines(new_ax)

            axes.append(new_ax)
        return axes


def uses_loaders(*loaders):
    """
    Class decorator to specify which data loaders a plot component requires.
    Updates the __required_loaders__ and __default_loader_kwargs__ attributes of the class.
    Also updates the __init__ method to include the loader arguments in the signature.
    If uses_loaders has been called multiple times, the loader arguments and the __init__ signature
    are overwritten.
    """
    def decorator(cls: Type[PlotComponent]):
        cls.__required_loaders__ = loaders
        loader_kwargs = _collect_all_kwargs(*loaders)
        cls.__loader_kwargs_signature__ = {**cls.__loader_kwargs_signature__, **loader_kwargs}
        if cls.__original_init__ is None:
            cls.__original_init__ = cls.__init__
        cls.__init__ = _update_signature(
            cls.__original_init__,
            cls.__loader_kwargs_signature__
        )
        return cls
    return decorator




# Helpers funcs
def _collect_all_kwargs(*loaders):
    loader_kwargs = {}
    for loader in loaders[::-1]:
        if not issubclass(loader, PlotDataLoader):
            raise ValueError(f"Loader {loader} is not a subclass of DataLoader.")
        loader_kwargs.update(loader.get_fullargspec())
    return loader_kwargs


def _update_signature(original_init, loader_kwargs: dict):
    """
    A decorator to dynamically update the init signature of a class
    by gathering parameters from all loaders.
    """
    original_signature = inspect.signature(original_init)
    signature = remove_signature_parameters(
        original_signature,
        'kwargs'
    )
    params = [
        inspect.Parameter(arg, inspect.Parameter.KEYWORD_ONLY, default=value)
        for arg, value in loader_kwargs.items()
    ]
    signature = add_signature_parameters(
        signature,
        last=params
    )

    signature = add_signature_parameters(
        signature,
        last=inspect.Parameter('plotting_kwargs', inspect.Parameter.VAR_KEYWORD)
    )
    @wraps(original_init, new_sig=signature)
    def wrapped_init(self, *args, **kwargs):
        return original_init(self, *args, **kwargs)
    
    return wrapped_init
