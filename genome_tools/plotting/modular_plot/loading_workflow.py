from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm

from genome_tools.plotting.modular_plot.utils import DataBundle
from genome_tools.plotting.modular_plot.display import markdown_value


class LoadingWorkflow:
    """
    Loader calls that `IntervalPlotter.get_interval_data` runs.
    Created by `IntervalPlotter.get_loading_workflow` (nothing is loaded), loaded with `load`.
    In notebooks shown as a tree: intervals, then loader calls with their args
    and the components that use them; `plot` draws it.

    Attributes
    ----------
    tree : nx.DiGraph
        `nx.prefix_tree` of the component paths (see `IntervalPlotter.resolve_loaders_caching`).
    calls : dict
        path element -> interval, or (component, LoaderClass, loader_kwargs).
    paths : list
        Path of each component.
    component_names : list
    last_nodes : list
        Tree node where the data of each component is complete.
    used_by : dict
        Tree node -> names of the components that use it.
    """
    def __init__(self, tree, calls, paths, component_names):
        self.tree = tree
        self.calls = calls
        self.paths = paths
        self.component_names = component_names
        self.used_by = {}
        self.last_nodes = []
        for name, path in zip(component_names, paths):
            node = 0
            for element in path:
                node = next(child for child in tree.successors(node) if tree.nodes[child]['source'] == element)
                self.used_by.setdefault(node, []).append(name)
            self.last_nodes.append(node)

    @classmethod
    def from_tasks(cls, tasks, use_loaders_cache=True):
        """
        Find the loader calls shared between components. Nothing is loaded.
        Each component is a path: its interval, then its loader calls. Paths with a common
        beginning share it in the prefix tree, so each shared loader call runs once.
        Intervals and loader kwargs are compared as objects (by id); `calls` holds them,
        so the compared ids stay unique while the workflow exists.
        If use_loaders_cache is False, every component gets its own path.

        Parameters
        ----------
        tasks : list of (component, interval, [(LoaderClass, loader_kwargs), ...])
            Output of `IntervalPlotter.resolve_loader_kwargs`.
        use_loaders_cache : bool, optional
        """
        calls = {}
        paths = []
        for i, (component, component_interval, resolved) in enumerate(tasks):
            start = ('interval', id(component_interval) if use_loaders_cache else i)
            calls[start] = component_interval
            path = [start]
            for LoaderClass, loader_kwargs in resolved:
                call = (LoaderClass, tuple((k, id(v)) for k, v in sorted(loader_kwargs.items())))
                calls.setdefault(call, (component, LoaderClass, loader_kwargs))
                path.append(call)
            paths.append(path)
        return cls(nx.prefix_tree(paths), calls, paths, [component.name for component, _, _ in tasks])

    def load(self, n_cpus=1):
        """
        Run the loader calls.
        With n_cpus > 1, each first loader call (with the calls that follow it) runs in its own process.

        Returns
        -------
        list of DataBundle
            Data of each component, in the order of self.component_names.
        """
        data = {}
        jobs = []
        for interval_node in self.tree.successors(0):
            data[interval_node] = DataBundle(interval=self.calls[self.tree.nodes[interval_node]['source']])
            # a job gets only its part of the tree and its calls (sent to another process with n_cpus > 1)
            for first in self.tree.successors(interval_node):
                if first == -1:  # end of a path
                    continue
                subtree = self.tree.subgraph({first, *nx.descendants(self.tree, first)}).copy()
                subtree_calls = {
                    subtree.nodes[node]['source']: self.calls[subtree.nodes[node]['source']]
                    for node in subtree if node != -1
                }
                jobs.append((subtree, subtree_calls, first, data[interval_node]))

        n_cpus = min(n_cpus, len(jobs))
        if n_cpus <= 1:
            for job in tqdm(jobs, desc="Loading data"):
                data.update(self._load_subtree(*job))
        else:
            with ProcessPoolExecutor(max_workers=n_cpus) as executor:
                futures = [executor.submit(self._load_subtree, *job) for job in jobs]
                for future in tqdm(as_completed(futures), total=len(futures), desc="Loading data"):
                    data.update(future.result())

        results = []
        handed_out = set()
        for node in self.last_nodes:
            # components with the same path get their own copy
            results.append(data[node].copy() if node in handed_out else data[node])
            handed_out.add(node)
        return results

    @staticmethod
    def _load_subtree(tree, calls, first, interval_data):
        """
        Run the loader call `first` on interval_data, then the calls that follow it.
        Returns {tree node: DataBundle}.
        """
        data = {}
        for node in nx.dfs_preorder_nodes(tree, first):
            if node == -1:  # end of a path
                continue
            component, LoaderClass, loader_kwargs = calls[tree.nodes[node]['source']]
            previous = interval_data if node == first else data[next(tree.predecessors(node))]
            # copy, so the result of the previous loader call is not modified
            data[node] = component.run_loader(previous.copy(), LoaderClass, loader_kwargs)
        return data

    @property
    def n_loader_calls(self):
        """
        Number of loader calls that run (shared calls count once).
        """
        intervals = set(self.tree.successors(0))
        return sum(1 for node in self.tree if node not in (0, -1) and node not in intervals)

    def __repr__(self):
        return f"LoadingWorkflow({len(self.paths)} components, {self.n_loader_calls} loader calls)"

    def _repr_markdown_(self):
        n_without_cache = sum(len(path) - 1 for path in self.paths)
        lines = [
            f"**LoadingWorkflow** · {len(self.paths)} components, {self.n_loader_calls} loader calls "
            f"({n_without_cache} without cache)",
            '',
        ]
        # one block per first loader call: its interval, then the calls that follow it,
        # so components that share calls are shown together; blocks in the order of the components
        blocks = []
        for interval_node in self.tree.successors(0):
            interval = self.calls[self.tree.nodes[interval_node]['source']]
            interval = interval.to_ucsc(commas=True) if hasattr(interval, 'to_ucsc') else repr(interval)
            depths = nx.shortest_path_length(self.tree, interval_node)
            for first in self.tree.successors(interval_node):
                if first == -1:  # components without loaders
                    names = [name for name, node in zip(self.component_names, self.last_nodes) if node == interval_node]
                    blocks.append((names, [f"- `{interval}` · {'; '.join(f'**{name}**' for name in names)}"]))
                    continue
                names = self.used_by[first]
                block = [f"- `{interval}` · {'; '.join(f'**{name}**' for name in names)}"]
                for node in nx.dfs_preorder_nodes(self.tree, first):
                    if node == -1:  # end of a path
                        continue
                    _, LoaderClass, loader_kwargs = self.calls[self.tree.nodes[node]['source']]
                    node_names = '; '.join(f"**{name}**" for name in self.used_by[node])
                    args = ', '.join(f"`{k}`={markdown_value(v)}" for k, v in loader_kwargs.items())
                    args = f" · {args}" if args else ''
                    # escaped number: a markdown ordered list would renumber the steps
                    block.append(f"  - {depths[node]}\\. `{LoaderClass.__name__}` · {node_names}{args}")
                blocks.append((names, block))
        for names, block in sorted(blocks, key=lambda b: self.component_names.index(b[0][0])):
            lines.extend(block)
        return '\n'.join(lines)

    def _plot_layout(self):
        """
        Graph to plot (prefix tree without its root and end-of-path markers),
        node positions and labels (interval or loader name, and the components that use it).
        Top-down tree: leaves side by side, parents centered above their children.
        """
        graph = self.tree.copy()
        graph.remove_nodes_from([0, -1])  # root and end-of-path markers of nx.prefix_tree
        intervals = list(self.tree.successors(0))
        labels = {}
        for node in graph:
            call = self.calls[graph.nodes[node]['source']]
            if node in intervals:
                name = call.to_ucsc(commas=True) if hasattr(call, 'to_ucsc') else repr(call)
            else:
                name = call[1].__name__
            labels[node] = name + '\nused by:\n' + '; '.join(self.used_by[node])

        # children are visited before their parent
        pos = {}
        n_leaves = 0
        for interval_node in intervals:
            depths = nx.shortest_path_length(graph, interval_node)
            for node in nx.dfs_postorder_nodes(graph, interval_node):
                children = list(graph.successors(node))
                if children:
                    x = sum(pos[child][0] for child in children) / len(children)
                else:
                    x = n_leaves
                    n_leaves += 1
                pos[node] = (x, -depths[node])
        return graph, pos, labels

    def plot(self, ax=None, edges_kwargs=None, **labels_kwargs):
        """
        Plot the workflow as a top-down tree: intervals at the top, then loader calls,
        each with the components that use it.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            If None, creates a figure sized to the tree.
        edges_kwargs : dict, optional
            Passed to `nx.draw_networkx_edges`, e.g. edge_color, width.
        **labels_kwargs
            Passed to `nx.draw_networkx_labels`, e.g. font_size (6), font_color, bbox.

        Returns
        -------
        ax : matplotlib.axes.Axes
        """
        graph, pos, labels = self._plot_layout()
        if ax is None:
            n_columns = 1 + max((x for x, _ in pos.values()), default=0)
            n_levels = 1 - min((y for _, y in pos.values()), default=0)
            _, ax = plt.subplots(figsize=(2.5 * n_columns, 1.2 * n_levels))
        edges_kwargs = {'arrows': False, **(edges_kwargs or {})}
        labels_kwargs = {
            'font_size': 6,
            'bbox': dict(boxstyle='round', facecolor='white', edgecolor='gray'),
            **labels_kwargs,
        }
        nx.draw_networkx_edges(graph, pos, ax=ax, **edges_kwargs)
        texts = nx.draw_networkx_labels(graph, pos, labels, ax=ax, **labels_kwargs)
        for text in texts.values():
            text.set_clip_on(False)
        ax.axis('off')
        return ax
