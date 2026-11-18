"""
Helpers for notebook display (`_repr_markdown_`) of components, data and loading workflows.
"""
import html


def markdown_value(value, max_len=50):
    """
    Short markdown summary of a value for notebook tables.
    Shapes instead of contents for arrays/DataFrames, lengths for long containers.
    A class can define `_markdown_value_()` to give its own summary (e.g. RequiredArgument).
    """
    if hasattr(type(value), '_markdown_value_'):
        return value._markdown_value_()
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


def markdown_parts(value):
    """
    One-line summary (inline HTML) and markdown body (None if nothing to expand) of a value.
    """
    if hasattr(value, '_markdown_parts'):
        return value._markdown_parts()
    return f"<code>{html.escape(repr(value))}</code>", None


def components_tuple_markdown(components_tuple):
    """
    Markdown of a namedtuple keyed by component names (see `utils.components_namedtuple`):
    each value as its name and one-line summary, collapsed if it has more to show
    (expanded if the tuple has a single value).
    """
    fields = components_tuple._fields
    example = fields[0] if fields else 'name'
    parts = [
        f"**{type(components_tuple).__name__}** · namedtuple with {len(components_tuple)} fields, "
        f"access as `.{example}`, `[0]` or `.get('{example}')`"
    ]
    for name, value in zip(fields, components_tuple):
        summary, body = markdown_parts(value)
        title = f"<b>{html.escape(name)}</b>: {summary}"
        if body is None:
            parts.append(title)
        elif len(components_tuple) == 1:
            parts.append(f"{title}\n\n{body}")
        else:
            parts.append(f"<details><summary>{title}</summary>\n\n{body}\n\n</details>")
    return '\n\n'.join(parts)
