import html
import logging
from collections import namedtuple


def markdown_value(value, max_len=50):
    """
    Short markdown summary of a value for notebook tables.
    Shapes instead of contents for arrays/DataFrames, lengths for long containers.
    """
    if isinstance(value, RequiredArgument):
        return '**required, not set**'
    shape = getattr(value, 'shape', None)
    if isinstance(shape, tuple):
        dtype = getattr(value, 'dtype', None)
        return f"{type(value).__name__} {' × '.join(map(str, shape))}" + (f" {dtype}" if dtype is not None else "")
    if isinstance(value, (list, tuple, dict, set)) and len(value) > 3:
        return f"{type(value).__name__} ({len(value)} items)"
    text = repr(value)
    if len(text) > max_len:
        half = (max_len - 1) // 2
        text = text[:half] + '…' + text[-half:]
    return '`' + text.replace('`', "'").replace('|', '\\|') + '`'


def markdown_table(header, rows):
    """
    Markdown table from a header and rows of already formatted cells.
    """
    lines = ['| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
    lines += ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows]
    return '\n'.join(lines)


def components_namedtuple(field_names):
    """
    Namedtuple class keyed by component names (data, axes, plot components).
    In notebooks each value is shown as its name and one-line summary, collapsed
    if it has more to show (expanded if the tuple has a single value).
    """
    base = namedtuple('ComponentNamesTuple', field_names)

    def get(self, name, default=None):
        """
        Value by component name, like dict.get.
        """
        return getattr(self, name) if name in self._fields else default

    def _repr_markdown_(self):
        example = self._fields[0] if self._fields else 'name'
        parts = [
            f"**{type(self).__name__}** · namedtuple with {len(self)} fields, "
            f"access as `.{example}`, `[0]` or `.get('{example}')`"
        ]
        for name, value in zip(self._fields, self):
            summary, body = _markdown_parts(value)
            title = f"<b>{html.escape(name)}</b>: {summary}"
            if body is None:
                parts.append(title)
            elif len(self) == 1:
                parts.append(f"{title}\n\n{body}")
            else:
                parts.append(f"<details><summary>{title}</summary>\n\n{body}\n\n</details>")
        return '\n\n'.join(parts)

    return type('ComponentNamesTuple', (base,), {'__slots__': (), 'get': get, '_repr_markdown_': _repr_markdown_})


def _markdown_parts(value):
    """
    One-line summary (inline HTML) and markdown body (None if nothing to expand) of a value.
    """
    if hasattr(value, '_markdown_parts'):
        return value._markdown_parts()
    return f"<code>{html.escape(repr(value))}</code>", None


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
