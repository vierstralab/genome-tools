from matplotlib import gridspec
from genome_tools.plotting.utils import clear_spines
from genome_tools.plotting import signal_plot

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders

from genome_tools.plotting.modular_plot.loaders.dhs import ComponentTracksLoader, DHSIndexLoader, DHSLoadingsLoader

from .abstract import SegmentPlotComponent


@uses_loaders(DHSIndexLoader)
class DHSIndexComponent(SegmentPlotComponent):
    """Plot DHS index segments within the interval.

    Loaders: ``DHSIndexLoader`` (inherits ``SegmentsLoader``)

    Required loader args:
    - ``dhs_index``: DataFrame of DHS rows overlapping the interval
      (extra columns may be attached as rectprops)

    Plot kwargs: forwarded to ``segment_plot``.

    Returns: ``matplotlib.axes.Axes``
    """
    __intervals_attr__ = DHSIndexLoader.__intervals_attr__

    def _plot(self, data, ax, **kwargs):
        super()._plot(data, ax, **kwargs)
        ax.set_ylabel('NMF-annotated\nDHSs')
        return ax
    

@uses_loaders(ComponentTracksLoader)
class NMFTracksComponent(IntervalPlotComponent):
    """Stacked DNase density tracks for NMF components.

    Loaders: ``ComponentTracksLoader``

    Required loader args:
    - ``cutcounts_files``: dict[int, list[str]] mapping component index -> bigWig files
    - ``nmf_components``: list of component indices to plot
    - optional ``smooth``, ``step``, ``bandwidth`` for aggregation

    Plot kwargs:
    - ``component_data``: DataFrame with columns ``index``, ``color``, ``name``
    - ``common_lim``: share y-limits across tracks

    Returns: tuple (top axes, list of density axes)
    """
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, component_data, **kwargs):
        density_axes = self.plot_component_tracks(data.interval, data.nmf_components,
                                                  data.component_tracks, component_data,
                                                  gridspec_ax=ax, **kwargs)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel("DNase I\ndensity")
        clear_spines(ax)
        return ax, density_axes
    

    @staticmethod
    def plot_component_tracks(interval, components, component_tracks, component_data,
                              gridspec_ax, common_lim=False, **kwargs):
        """Helper to render per-component density subplots.

        Parameters:
        - interval: GenomicInterval
        - components: list[int]
        - component_tracks: list[np.ndarray] densities per component
        - component_data: pd.DataFrame with color/name per component
        - gridspec_ax: parent Axes to host subplots
        - common_lim: bool to unify y-limits
        - **kwargs: forwarded to ``signal_plot``

        Returns: list[matplotlib.axes.Axes]
        """
        assert len(components) == len(component_tracks)
        gss = gridspec.GridSpecFromSubplotSpec(len(components), 1, subplot_spec=gridspec_ax, hspace=0.05)
        axes = []
        colors = component_data.set_index('index').loc[components, 'color'].values
        labels = component_data.set_index('index').loc[components, 'name'].values
        lims = []
        for i, (color, label, segs) in enumerate(zip(colors, labels, component_tracks)):
            ax_dens = gridspec_ax.get_figure().add_subplot(gss[i])
            lim = segs.max()
            lims.append(lim)
            signal_plot(interval, segs, ax=ax_dens, color=color, lw=0, **kwargs)
            ax_dens.set_ylim(0, lim)
            ax_dens.set_yticks([])
            if i != len(components) - 1:
                ax_dens.set_xticks([])
            ax_dens.text(x=0.005, y=0.45, ha='left', va='bottom', s=label, transform=ax_dens.transAxes, fontsize=5)
            axes.append(ax_dens)

        if common_lim:
            for ax in axes:
                ax.set_ylim(0, max(lims))
        return axes
    

@uses_loaders(DHSIndexLoader, DHSLoadingsLoader)
class DHSLoadingsComponent(IntervalPlotComponent):
    """Per-DHS NMF loadings barplots centered at DHS summits.

    Loaders: ``DHSIndexLoader``, ``DHSLoadingsLoader``

    Required loader args:
    - ``dhs_index``: DataFrame of DHS rows
    - ``H``: NMF loadings matrix (components x DHS)

    Plot kwargs:
    - ``component_data``: DataFrame describing components (color/name)
    - ``bp_width``: width around each summit to display

    Returns: tuple (top axes, list of per-DHS axes)
    """

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, component_data, bp_width=50, **kwargs):
        ax.axis('off')
        dhs_intervals = getattr(data, DHSIndexLoader.__intervals_attr__)
        axes = self.add_axes_at_summits(
            dhs_intervals,
            data.interval,
            ax=ax,
            bp_width=bp_width
        )
        self.plot_barplots_for_dhs(dhs_intervals, axes, H=data.H, component_data=component_data)
        return ax, axes
    
    @staticmethod
    def plot_barplots_for_dhs(genomic_intervals, axes, H, component_data):
        assert len(genomic_intervals) == len(axes)
        from nmf_tools.plotting.matrices_barplots import component_barplot

        for genomic_interval, ax in zip(genomic_intervals, axes):
            component_barplot(H[:, genomic_interval.index: genomic_interval.index + 1], component_data, ax=ax, normalize=True)


