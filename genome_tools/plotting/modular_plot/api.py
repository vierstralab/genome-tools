import inspect
import sys

from matplotlib import pyplot as plt
from typing import List, Type, Union, Sequence
import types
from makefun import wraps, add_signature_parameters, remove_signature_parameters, with_signature
import warnings

from genome_tools.plotting.utils import format_axes_to_interval, add_axes_at_intervals
from genome_tools import GenomicInterval

from genome_tools.plotting.modular_plot.utils import LoggerMixin, DataBundle, RequiredArgument


class PlotDataLoader(LoggerMixin):
    """
    Base class for all data loaders.
    Each loader should implement the `_load` method to load and filter data based on the interval.

    By default (if not implemented), the _load method is created from the required_loader_kwargs. The generated _load method takes required_loader_kwargs as arguments and sets them as fields of the data object.
    """
    required_loader_kwargs = []
    uses_custom_load = True

    def __init__(self, logger_level=None):
        self.name = self.__class__.__name__
        super().__init__(logger_level=logger_level)
    
    def __repr__(self):
        if self.uses_custom_load:
            return f"{self.name} with custom _load method"
        return f"{self.name} sets fields: {', '.join(self.required_loader_kwargs)}"

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls._load is PlotDataLoader._load:
            cls._set_default_load()
            cls.uses_custom_load = False
        elif len(cls.required_loader_kwargs) > 0:
            warnings.warn(
                f"Both required_loader_kwargs and _load are specified for loader {cls.__name__}. required_loader_kwargs argument is ignored.",
                stacklevel=3
            )

    @classmethod
    def _set_default_load(cls):
        """
        Set the default load method if not already set.
        """
        if len(cls.required_loader_kwargs) == 0:
            warnings.warn(
                f"Loader {cls.__name__} has no required_loader_kwargs and no _load method implemented. The loader does not modify the data.",
                stacklevel=3
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
    __loader_arg_sources__ = {} # maps loader args to the loaders that use them
    __original_init__ = None

    def __init__(self, name=None, logger_level=None, **kwargs):
        super().__init__(logger_level=logger_level)
        if name is None:
            name = self.__class__.__name__
        self.name = name

        # Separate loader kwargs from plot kwargs
        self.loader_overrides = self._extract_loader_specific_overrides(
            kwargs,
            loaders=self.__required_loaders__
        )
        self.loader_kwargs = {
            k: kwargs.pop(k) for k in list(kwargs)
            if k in self.__class__.__loader_kwargs_signature__
        }
        self._validate_loader_arg_sources()

        self.plot_kwargs = kwargs

    def __repr__(self):
        return f"{self.__class__.__name__}(name={self.name}, loaders={[loader.__name__ for loader in self.__required_loaders__]})"
    
    def _validate_loader_arg_sources(self):
        """
        Validate that the loader arg sources are consistent with the required loaders.
        """
        for arg, entry in self.__loader_arg_sources__.items():
            if len(entry['loaders']) > 1:
                overlapping_loaders = ', '.join([loader.__name__ for loader in entry['loaders']])
                val = {arg: 'some_value'}
                example = f"{entry['loaders'][0].__name__}={val}"
                self.logger.warning(
                    f"Argument '{arg}' is used by multiple loaders "
                    f"({overlapping_loaders}). "
                    f"The same value will be passed to all of them.\nOr you can specify loader-specific overrides as follows ({example})."
                )
    
    @staticmethod
    def _extract_loader_specific_overrides(kwargs: dict, classes: List[Union[PlotDataLoader, 'PlotComponent']]):
        """
        Extract loader or component specific override dicts from kwargs.

        Parameters
        ----------
        kwargs : dict
            The keyword arguments passed to the plot component or loader.
        
        classes : List[Union[PlotDataLoader, 'PlotComponent']]
            The list of loader or component classes to extract overrides for.

        Example:
            kwargs = {"SomeLoader": {"a": 1}, "color": "red"}
            returns {"SomeLoader": {"a": 1}}
            and kwargs becomes {"color": "red"}

        Returns
        -------
        dict : Dict[str, dict]
            A dictionary mapping loader class names to their specific override dictionaries.
        """
        overrides = {}

        for key in list(kwargs.keys()):
            for class_ in classes:
                name = class_.name
                if key == name:
                    val = kwargs.pop(key)
                    if not isinstance(val, dict):
                        raise TypeError(
                            f"Loader-specific override for '{key}' must be a dict, "
                            f"got {val} {type(val).__name__}"
                        )
                    overrides[name] = val
                    break

        return overrides

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
        new_class = uses_loaders(*loaders)(new_class)
        module = sys.modules[cls.__module__]

        # Protect against name collisions
        if hasattr(module, new_class_name):
            raise RuntimeError(
                f"Class name '{new_class_name}' already exists in module {cls.__module__}."
                " Choose a different new_class_name."
            )

        # Register class so pickle/multiprocessing can find it
        setattr(module, new_class_name, new_class)
        new_class.__module__ = cls.__module__
        return new_class
    
    def load_data(self, data: DataBundle, **runtime_overrides):
        """
        Modifies data for the plot component using the required loaders.

        Merge priority:
            1. loader defaults     (lowest)
            2. component overrides (self.loader_kwargs)
            3. runtime overrides   (**runtime_loader_kwargs)
            4. component loader-specific overrides (self.loader_overrides)
            5. runtime loader-specific overrides (_extract_loader_specific_overrides(runtime_overrides)) (highest)
        """
        for LoaderClass in self.__required_loaders__:
            loader_defaults = LoaderClass.get_fullargspec()
            component_overrides = self.loader_kwargs
            component_loader_specific_overrides = self.loader_overrides.get(LoaderClass.__name__, {})
            runtime_loader_specific_overrides = self._extract_loader_specific_overrides(
                runtime_overrides,
                loaders=[LoaderClass]
            )

            kwargs_for_loader = {
                **loader_defaults,
                **component_overrides, # specify during component init
                **runtime_overrides, # specify during load_data_for_interval
                **component_loader_specific_overrides, # specify during component init
                **runtime_loader_specific_overrides, # specify during load_data_for_interval
            }

            # Keep only args that loader actually accepts
            kwargs_for_loader = {
                k: v for k, v in kwargs_for_loader.items()
                if k in loader_defaults
            }

            loader = LoaderClass(
                logger_level=self.logger.level
            )
            data = loader.load(data, **kwargs_for_loader)
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
                 name=None,
                 height: float = 1.0, 
                 margins: Union[float, Sequence[float]] = 0.1,
                 interval_key=None,
                 logger_level=None,
                 **kwargs):
        self.height = height
        self.margin_top, self.margin_bottom = self._parse_margins(margins)
        self.interval_key = interval_key

        super().__init__(name=name, logger_level=logger_level, **kwargs)

    @property
    def full_height(self):
        """Calculate the full height of the plot including margins."""
        return self.height + self.margin_top + self.margin_bottom

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
            format_axes_to_interval(ax, data.interval, axis='x')
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
        return add_axes_at_intervals(
            modified_genomic_intervals,
            interval,
            ax=ax
        )


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
        loader_kwargs, args_info = _collect_all_kwargs(*loaders)
        cls.__loader_arg_sources__ = args_info
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
    args_info = {}
    for loader in loaders[::-1]:
        if not issubclass(loader, PlotDataLoader):
            raise ValueError(f"Loader {loader} is not a subclass of PlotDataLoader.")
        
        # Default is not really used
        spec = loader.get_fullargspec()
        for arg, default in spec.items():
            entry = args_info.get(
                arg,
                {"default": default, "loaders": []}
            )
            if isinstance(entry["default"], RequiredArgument):
                entry["default"] = default
            entry["loaders"].append(loader)
            args_info[arg] = entry

    loader_kwargs = {arg: entry["default"] for arg, entry in args_info.items()}

    return loader_kwargs, args_info


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
