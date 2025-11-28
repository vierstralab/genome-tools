from nbformat import reads
import numpy as np
import pandas as pd
from typing import List

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.offsetbox import AnnotationBbox, TextArea, HPacker

from genome_tools import VariantInterval, GenomicInterval
from genome_tools.plotting.colors.cm import get_vocab_color
from genome_tools.plotting import segment_plot
from genome_tools.plotting.sequence import plot_letter
from genome_tools.plotting.utils import format_axes_to_interval, pack_rows, add_axes_at_intervals

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.loaders.variant import FinemapLoader, AggregatedCAVLoader, PerSampleCAVLoader, AllelicReadsLoaderFPTools, VariantGenotypeLoader


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
    def plot_single_bp_objects(variant_intervals, interval, ax=None, s=1.0, **kwargs):
        if ax is None:
            ax = plt.gca()

        positions = [v.pos - 0.5 for v in variant_intervals]
        values = [v.value for v in variant_intervals]
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
    def _plot(self, data, ax, **kwargs):

        self.plot_single_bp_objects(
            data.cavs_intervals,
            data.interval,
            ax=ax
        )
        self.annotate_variant_alleles(data.cavs_intervals, ax=ax, **kwargs)
        max_value = max(data.cavs_intervals, key=lambda x: x.value).value
        ax.set_ylim(0, max(max_value * 1.2, 2.0))
        ax.set_xticks([])
        ax.set_ylabel("CAV\neffect size")
        return ax
    

NonAggregatedCAVComponent = CAVComponent.with_loaders(PerSampleCAVLoader, new_class_name="NonAggregatedCAVComponent")


@uses_loaders(VariantGenotypeLoader, AllelicReadsLoaderFPTools)
class AllelicCutcountsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax: plt.Axes, **kwargs):

        ref_cuts = data.ref_cutcounts
        alt_cuts = data.alt_cutcounts
        interval: VariantInterval = data.variant_interval

        tot_reads = ref_cuts.sum() + alt_cuts.sum()
        
        # Shared x-axis

        if ref_cuts.sum() > alt_cuts.sum():
            ylim_ref = np.quantile(ref_cuts, 0.75)
            ylim_alt = ylim_ref/3 #np.quantile(alt_y, 0.65)
        else:
            ylim_alt = np.quantile(alt_cuts, 0.75)
            ylim_ref = ylim_alt/3 #np.quantile(ref_y, 0.65)

        # Plotting
        gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=ax.get_subplotspec(), hspace = 0.25)
        ax.axis("off")
        fig = ax.get_figure()
        # fig.suptitle(data.rs_id)
        axes = []
        x = np.arange(data.interval.start, data.interval.end) + 0.5 # centering bars
        for i, (cuts, allele, ylim) in enumerate(zip(
            [ref_cuts, alt_cuts],
            [interval.ref, interval.alt],
            [ylim_ref, ylim_alt]
        )):
            ax_bar = fig.add_subplot(gs[i, :])
            ax_bar.bar(
                x,
                cuts,
                width=1,
                color=get_vocab_color(allele, 'dna'), 
                label=f"{allele}: {cuts.sum()}\n({round(cuts.sum() / tot_reads * 100, 2)}%)"
            )
            ax_bar.set_ylim(0, ylim)
            ax_bar.xaxis.set_visible(False)
            format_axes_to_interval(ax_bar, data.interval, axis='x')
            axes.append(ax_bar)
        ax.legend(loc='upper right', fontsize='small')
                    
        return ax, axes


@uses_loaders(VariantGenotypeLoader, AllelicReadsLoaderFPTools)
class AllelicReadsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, reads_count_tr=120, max_distance_to_read_ends=200, pad_bp=5, rect_height=0.4, **kwargs):
        ref_reads: List[GenomicInterval] = data.ref_reads#, key=lambda x: x.start + x.end)
        alt_reads: List[GenomicInterval] = data.alt_reads#, key=lambda x: x.start + x.end)
        variant_interval: VariantInterval = data.variant_interval
        for r in ref_reads:
            r.rectprops = dict(color=get_vocab_color(r.base, 'dna', default='grey'))
        for r in alt_reads:
            r.rectprops = dict(color=get_vocab_color(r.base, 'dna', default='grey'))
        reads = ref_reads + alt_reads


        reads = [
            read for read in reads 
            if (
                abs(read.start - variant_interval.start) <= max_distance_to_read_ends and
                abs(read.end - variant_interval.end) <= max_distance_to_read_ends
            )
        ]

        if len(reads) > reads_count_tr: # subsample
            rng = np.random.default_rng(seed=42)
            idx = rng.choice(len(reads), size=reads_count_tr, replace=False)
            idx.sort()
            reads = [reads[i] for i in idx]

        _, reads = pack_rows(reads, pad_bp=pad_bp)
        letter_pad = (1.0 - rect_height) / 2.0 / 2.0
        letter_height = 1 - 2 * letter_pad
        letter_width = 0.7

        # Add axes for letters
        letter_axes = add_axes_at_intervals(
            reads,
            data.interval,
            ax=ax,
        )
        for read, letter_ax in zip(reads, letter_axes):
            plot_letter(
                letter=read.base,
                x=1.0 - letter_width / 2.0,
                y=letter_pad,
                height=letter_height,
                width=letter_width,
                ax=letter_ax,
            )
            letter_ax.set_xlim(0, 1.0)
            letter_ax.set_ylim(0, 1.0)
            letter_ax.patch.set_color('white')
            letter_ax.set_xticks([])
            letter_ax.set_yticks([])
            for s in letter_ax.spines.values():
                s.set_visible(True)
                s.set_linewidth(0.5)
            

        segment_plot(data.interval, reads, rect_height=rect_height, pack=False, ax=ax, **kwargs)
        
        ax.set_yticks([])
        return ax