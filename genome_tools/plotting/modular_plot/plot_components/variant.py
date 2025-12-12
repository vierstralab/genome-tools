import numpy as np
from typing import List

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.offsetbox import AnnotationBbox, TextArea, HPacker
import matplotlib.patches as mpatches

from genome_tools import VariantInterval, GenomicInterval
from genome_tools.plotting.colors.cm import get_vocab_color
from genome_tools.plotting import segment_plot
from genome_tools.plotting.sequence import plot_letter
from genome_tools.plotting.utils import format_axes_to_interval, pack_rows, add_axes_at_intervals

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.loaders.variant import (
    FinemapLoader,
    AggregatedCAVLoader,
    PerSampleCAVLoader,
    ReadsLoader,
    AllelicReadsLoader,
    VariantGenotypeLoader,
    VariantIntervalLoader,
    GroupsByGenotypeLoader
)
from genome_tools.plotting.modular_plot.loaders.footprint import GroupsDataLoader


# TODO fix other components
class LolipopVariantsComponent(IntervalPlotComponent):
    """
    A vertical plot component that plots single base pair objects
    within a genomic interval.
    self.plot accepts a data object with 'positions' and 'values' fields
    to plot them on the x-axis and y-axis, respectively as scatter points
    and vertical lines from the x-axis.

    kwargs are passed to the scatter plot function.
    """
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        self.plot_single_bp_objects(data.positions, data.values, data.interval, ax=ax, **kwargs)
        return ax
    
    @staticmethod
    def plot_single_bp_objects(variant_intervals, interval, abs_es=True, ax=None, s=1.0, **kwargs):
        if ax is None:
            ax = plt.gca()

        positions = [v.pos - 0.5 for v in variant_intervals]
        values = [v.value for v in variant_intervals]
        if abs_es:
            values = [abs(v) for v in values]
        colors = [getattr(v, 'color', 'k') for v in variant_intervals]
        annotations = [getattr(v, 'annotation', None) for v in variant_intervals]

        ax.scatter(x=positions, y=values, s=s, c=colors, **kwargs)
        for val, pos, annotation, color in zip(values, positions, annotations, colors):
            ax.plot([pos, pos], [0, val], lw=0.5, color=color)
            if annotation is not None:
                ax.text(s=annotation, x=pos + (interval.end - interval.start) * 0.01,
                        y=val, fontsize=5, ha='left', va='center')
        if all(val >= 0 for val in values):
            ax.set_ylim(0, max(values) * 1.1)
        return ax
    
    @staticmethod
    def annotate_variant_alleles(variant_intervals, ax=None, box_alignment=(0.5, -0.6), **kwargs):
        if ax is None:
            ax = plt.gca()

        for v in variant_intervals:
            x = v.pos - 0.5
            y = v.value
            a1, a2 = v.ref, v.alt

            text_areas = []
            for text, color in zip([a1, '/', a2], [get_vocab_color(a1, 'dna'), 'k', get_vocab_color(a2, 'dna')]):
                text_area = TextArea(text, textprops=dict(color=color, fontsize=5, ha='center', va='center'))
                text_areas.append(text_area)
            
            hp = HPacker(children=text_areas, align="bottom", pad=0, sep=1)
            ab = AnnotationBbox(hp, (x, y), frameon=False, box_alignment=box_alignment, xycoords='data', **kwargs)
            ax.add_artist(ab)
            
        return ax



@uses_loaders(FinemapLoader)
class FinemapComponent(LolipopVariantsComponent):
    def _plot(self, data, ax, **kwargs):
        data.positions = data.unique_finemap_df['start'] + 0.5
        data.values = data.unique_finemap_df['pip']
        return super().plot(data, ax, **kwargs)


@uses_loaders(AggregatedCAVLoader)
class CAVComponent(LolipopVariantsComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, abs_es=True, **kwargs):

        self.plot_single_bp_objects(
            data.cavs_intervals,
            data.interval,
            abs_es=abs_es,
            ax=ax
        )
        self.annotate_variant_alleles(data.cavs_intervals, ax=ax, **kwargs)
        max_value = max(data.cavs_intervals, key=lambda x: x.value).value
        ax.set_ylim(0, max(max_value * 1.2, 2.0))
        ax.set_xticks([])
        ax.set_ylabel("CAV\neffect size")
        return ax
    

NonAggregatedCAVComponent = CAVComponent.with_loaders(PerSampleCAVLoader, new_class_name="NonAggregatedCAVComponent")

@uses_loaders(GroupsDataLoader, ReadsLoader)
class CutcountsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax: plt.Axes, vocab='dna', window=1, **kwargs):
        from hotspot3.connectors.bottleneck import BottleneckWrapper

        tracks = data.cutcount_tracks
        bn = BottleneckWrapper()

        # Total reads for labeling
        tot_reads = sum(t["cutcounts"].sum() for t in tracks)

        # Plotting
        gs = gridspec.GridSpecFromSubplotSpec(
            len(tracks),
            1,
            subplot_spec=ax.get_subplotspec(), 
            hspace=0.25
        )
        ax.axis("off")
        fig = ax.get_figure()

        axes = []
        x = np.arange(data.interval.start, data.interval.end) + 0.5

        for i, track in enumerate(tracks):
            
            cuts: np.ndarray = track["cutcounts"]
            aggregated_cuts = bn.centered_running_nanmean(cuts, window=window)
            allele: str = track["allele"]

            ax_bar = fig.add_subplot(gs[i, :])
            ax_bar.bar(
                x,
                aggregated_cuts,
                width=1,
                color=get_vocab_color(allele, vocab),
                label=f"{allele}: {cuts.sum()} ({round(cuts.sum() / tot_reads * 100, 2)}%)",
                **kwargs
            )
            ax_bar.xaxis.set_visible(False)
            format_axes_to_interval(ax_bar, data.interval, axis='x')
            axes.append(ax_bar)
        return ax, axes


AllelicCutcountsComponent = CutcountsComponent.with_loaders(VariantIntervalLoader, VariantGenotypeLoader, GroupsByGenotypeLoader, AllelicReadsLoader, new_class_name="AllelicCutcountsComponent")


@uses_loaders(VariantIntervalLoader, VariantGenotypeLoader, GroupsByGenotypeLoader, AllelicReadsLoader)
class AllelicReadsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, reads_count_tr=120, max_distance_to_read_ends=200, pad_bp=5, rect_height=0.4, circle_size=0.9, letter_size=0.5, font='IBM Plex Mono', seed=42, **kwargs):
        reads: List[GenomicInterval] = data.reads
        variant_interval: VariantInterval = data.variant_interval
        for r in reads:
            r.rectprops = dict(color=get_vocab_color(r.base, 'dna', default='grey'))

        reads = [
            read for read in reads 
            if (
                abs(read.start - variant_interval.start) <= max_distance_to_read_ends
                and abs(read.end - variant_interval.end) <= max_distance_to_read_ends 
                #and len(read) >= 10
            )
        ]

        if len(reads) > reads_count_tr: # subsample
            rng = np.random.default_rng(seed=seed)
            idx = rng.choice(len(reads), size=reads_count_tr, replace=False)
            idx.sort()
            reads = [reads[i] for i in idx]

        _, reads = pack_rows(reads, pad=pad_bp)

        letter_x_size = letter_size
        letter_y_size = letter_size

        nrows = max(read.row_index for read in reads) + 1
        dx_bp = self.frac_axis_height_to_bp_x(
            frac=1.0 / nrows,
            ax=ax,
        )

        letter_intervals = [
            GenomicInterval(
                chrom=variant_interval.chrom,
                start=variant_interval.start + 0.5 - dx_bp / 2.0,
                end=variant_interval.end + 0.5 + dx_bp / 2.0,
                row_index=read.row_index,
                base=read.base,
            ) for read in reads
        ]

        # Add axes for letters
        letter_axes = add_axes_at_intervals(
            letter_intervals,
            data.interval,
            row_pad=0.0,
            ax=ax,
        )
        for letter_interval, letter_ax in zip(letter_intervals, letter_axes):
            circ = mpatches.Circle(
                (0.5, 0.5),
                0.5 * circle_size,
                facecolor='white',
                edgecolor='k',
                linewidth=0.3,
                transform=letter_ax.transAxes,
                clip_on=False,
            )
            circ.set_clip_on(False)
            letter_ax.add_patch(circ)
            plot_letter(
                letter=letter_interval.base,
                x=0.0,
                y=(1 - letter_y_size) / 2,
                height=letter_y_size,
                width=letter_x_size,
                center_scale=True,
                font=font,
                ax=letter_ax,
            )
            letter_ax.set_xlim(0, 1.0)
            letter_ax.set_ylim(0, 1.0)
            letter_ax.set_xticks([])
            letter_ax.set_yticks([])
            letter_ax.axis('off')

        segment_plot(data.interval, reads, rect_height=rect_height, pack=False, ax=ax, **kwargs)
        
        ax.set_yticks([])
        ax.axis('off')
        return ax, letter_axes

    @staticmethod
    def frac_axis_height_to_bp_x(frac, ax):
        fig = ax.figure
        fig.canvas.draw_idle()

        bbox = ax.get_window_extent()
        height_px = bbox.height
        disp_len_px = frac * height_px

        x0, x1 = ax.get_xlim()
        x_ref = 0.5 * (x0 + x1)
        y0, y1 = ax.get_ylim()
        y_ref = 0.5 * (y0 + y1)

        trans = ax.transData
        inv   = trans.inverted()

        p0_disp = trans.transform((x_ref, y_ref))
        p1_disp = (p0_disp[0] + disp_len_px, p0_disp[1])

        (x1_data, _) = inv.transform(p1_disp)

        dx_bp = x1_data - x_ref
        return dx_bp
