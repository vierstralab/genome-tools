import html
import logging
from collections import namedtuple

from genome_tools.plotting.modular_plot.display import markdown_value, components_tuple_markdown


def components_namedtuple(field_names):
    """
    Namedtuple class keyed by component names (data, axes, plot components).
    In notebooks each value is shown as its name and one-line summary, collapsed
    if it has more to show (expanded if the tuple has a single value).
    """
    base = namedtuple('ComponentNamesTuple', field_names)
    return type('ComponentNamesTuple', (base,), {
        '__slots__': (),
        'get': _components_tuple_get,
        '_repr_markdown_': components_tuple_markdown,
    })


def _components_tuple_get(self, name, default=None):
    """
    Value by component name, like dict.get.
    """
    return getattr(self, name) if name in self._fields else default


class LoggerMixin:
    def __init__(self, logger_level=None):
        if logger_level is None:
            logger_level = logging.INFO
        self.logger = logging.getLogger(self.__class__.__name__)
        if not self.logger.hasHandlers():
            self.logger.setLevel(logger_level)

class DataBundle(LoggerMixin):
    """
    A container for holding all data needed by plot components.

    Attributes are added dynamically by the loaders
    depending on the data needed by the plot components.
    """
    def __init__(self, **kwargs):
        object.__setattr__(self, '_initialized_attributes', set())
        # Add additional data fields here (not strictly necessary)

        for key, value in kwargs.items():
            setattr(self, key, value)

        self.processed_loaders = []
        super().__init__()
    
    def __getattr__(self, name):
        loaders = ', '.join([loader.__name__ for loader in self.processed_loaders])
    
        raise AttributeError(f"Data loaded with {loaders} is missing attribute '{name}'. Available attributes: {self._get_display_attributes()}. Did you forget to attach the loader that sets this attribute?")
    
    def __setattr__(self, name, value):
        if name != 'logger' and name in self._initialized_attributes:
            self.logger.warning(f"Warning: Attribute '{name}' is overridden.")
        self._initialized_attributes.add(name)
        super().__setattr__(name, value)
    
    def __getitem__(self, key):
        return getattr(self, key)
    
    def _get_display_attributes(self):
        return {
            attr: getattr(self, attr) 
            for attr in self._initialized_attributes 
            if attr not in ('logger', 'processed_loaders', 'interval')
        }

    def __str__(self):
        keys = ', '.join(self._get_display_attributes().keys())
        return f"DataBundle({keys})"

    def __repr__(self):
        MAX_REPR_LEN = 50
        repr_str = f"DataBundle({self._get_display_attributes()})"
        if len(repr_str) > MAX_REPR_LEN:
            return str(self)
        return repr_str

    def _repr_markdown_(self):
        summary, body = self._markdown_parts()
        return summary if body is None else f"{summary}\n\n{body}"

    def _markdown_parts(self):
        """
        One-line summary (inline HTML) and markdown body (None if no fields) for notebook display.
        """
        summary = 'DataBundle'
        if 'interval' in self._initialized_attributes:
            interval = self.interval
            interval = interval.to_ucsc(commas=True) if hasattr(interval, 'to_ucsc') else repr(interval)
            summary += f" <code>{html.escape(interval)}</code>"
        fields = self._get_display_attributes()
        if not fields:
            return summary, None
        example = sorted(fields)[0]
        body = '\n'.join(f"- `{attr}`: {markdown_value(fields[attr])}" for attr in sorted(fields))
        body += f"\n\nAccess as `.{example}` or `['{example}']`"
        return summary, body

    def __getstate__(self):
        # Copy minimal state needed for multiprocessing
        state = self.__dict__.copy()
        # Avoid walking loader classes (causes recursion)
        state["processed_loaders"] = []
        return state
    
    def __setstate__(self, state):
        self.__dict__.update(state)

    def copy(self):
        """
        Create a shallow copy of the data bundle.
        """
        data = DataBundle()
        for attr in self._initialized_attributes:
            if attr in ('logger', 'processed_loaders'):
                continue
            setattr(data, attr, getattr(self, attr))
        data.processed_loaders.extend(self.processed_loaders)
        return data


class RequiredArgument:
    def __repr__(self):
        return 'Required loader arg'

    def _markdown_value_(self):
        return '**required, not set**'
