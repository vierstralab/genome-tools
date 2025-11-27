import numpy as np
import pandas as pd
import scipy
from scipy.special import logsumexp
import scipy.stats as st
import anndata as ad

from genome_tools.data.extractors import TabixExtractor
from genome_tools import GenomicInterval, VariantInterval, filter_df_to_interval


from footprint_tools.cli.post import posterior_stats as PosteriorStats
from footprint_tools.stats.windowing import stouffers_z
from footprint_tools.stats import posterior, differential
from footprint_tools.modeling import dispersion
from footprint_tools.stats.distributions import invchi2


from genome_tools.plotting.modular_plot import PlotDataLoader
from genome_tools.plotting.modular_plot.utils import DataBundle

from .basic import SegmentsLoader


class FootprintsIndexLoader(SegmentsLoader):
    __intervals_attr__ = 'footprint_intervals'
    def _load(self, data: DataBundle, footprints_index: pd.DataFrame, extra_columns=None, rectprops_columns=None):
        return super()._load(
            data, 
            segments_df=footprints_index, 
            extra_columns=extra_columns,
            rectprops_columns=rectprops_columns
        )


class PosteriorLoader(PlotDataLoader):
    def _load(
            self,
            data: DataBundle,
            posterior_file, 
            footprints_metadata: pd.DataFrame, 
            sort_heatmap_by_region: GenomicInterval = None,
            grouping_column='extended_annotation'
    ):
        # Get posterior data
        with TabixExtractor(
            posterior_file,
            skiprows=1
        ) as extractor:
            raw_posterior = extractor[data.interval].rename(columns={'# chrom': '#chr'})

        interval_posterior = raw_posterior.set_index('start').drop(
            columns=['#chr', 'end']
        ).T.loc[footprints_metadata.index]
        
        interval_posterior_df = pd.DataFrame(
            0.0,
            index=interval_posterior.index,
            columns=np.arange(data.interval.start, data.interval.end)
        )
        interval_posterior_df.loc[:, interval_posterior.columns] = interval_posterior
        
        grouping_column_data = footprints_metadata.loc[interval_posterior_df.index, grouping_column]
        if sort_heatmap_by_region is not None:
            order = self.sort_by_interval(
                interval_posterior_df,
                data.interval,
                sort_heatmap_by_region,
                grouping_column_data
            )

            interval_posterior_df = interval_posterior_df.loc[order]
            grouping_column_data = grouping_column_data.loc[order]
            
        data.interval_posterior = 1 - np.exp(-interval_posterior_df)
        data.grouping_column = grouping_column_data
        return data
    
    @staticmethod
    def sort_by_interval(df: pd.DataFrame, base_interval: GenomicInterval, region: GenomicInterval, group_column: pd.Series):
        assert region.chrom == base_interval.chrom
        assert region.start >= base_interval.start and region.end <= base_interval.end

        region_slice = slice(region.start - base_interval.start, region.end - base_interval.start)
        sample_means = df.loc[:, region_slice].mean(axis=1)
        group_mean = sample_means.groupby(group_column).transform("mean")
        order = pd.DataFrame({
            'group_mean': group_mean,
            'per_sample': sample_means,
            'group_col': group_column
        }, index=df.index).sort_values(['group_mean', 'per_sample'], ascending=True).index
        return order


class ProtectedNucleotidesLoader(PlotDataLoader):
    
    def _load(self, data: DataBundle, protected_nuc_sample_ids=None, threshold=0.99):
        """
        Process the posterior probability data
        """
        interval_posterior = data.interval_posterior # data.interval_posterior
        sequence = data.sequence
        
        if protected_nuc_sample_ids is not None:
            interval_posterior = interval_posterior.loc[protected_nuc_sample_ids, :]

        binary_df = (interval_posterior >= threshold).astype(int)
        col_sums = binary_df.sum(axis=0)
        # normalized based on max value
        normalized = col_sums / col_sums.max()
        # Sequence
        data.matrix = self.seq_heights_to_matrix(sequence, normalized)
        return data

    @staticmethod
    def seq_heights_to_matrix(seq, heights):
        """
        Convert a sequence string and a height vector (Series/array)
        into a (L x 4) matrix suitable for plotting a sequence logo.
        """
        vocab = ["A", "C", "G", "T"]
        letter_to_index = {l: i for i, l in enumerate(vocab)}
        mat = np.zeros((len(seq), len(vocab)))
        for i, (base, h) in enumerate(zip(seq, heights)):
            if base in letter_to_index:
                mat[i, letter_to_index[base]] = h
        return mat #mat.T  


class FootprintsDataLoader(PlotDataLoader):
    def _load(self, data: DataBundle, footprints_metadata: pd.DataFrame, variant_interval: VariantInterval, calc_posteriors=False):
        variant_genotypes: pd.DataFrame = data.variant_genotypes # indiv_id, variant pairs
        variant_genotypes = filter_df_to_interval(variant_genotypes, variant_interval)

        samples_with_genotype = footprints_metadata.dropna(
            subset='indiv_id'
        ).reset_index(
            names='sample_id'
        ).merge(
            variant_genotypes
        ).set_index(
            "sample_id"
        ).sort_values(
            by="parsed_genotype"
        )

        data.a_allele_count = sum(samples_with_genotype['parsed_genotype'] == 'A')

        data.disp_models = [
            dispersion.load_dispersion_model(x) 
            for x in samples_with_genotype["dm_file"].values
        ]

        obs, exp, fdr, w = self.load_fp_data(
            samples_with_genotype["tabix_file"].tolist(),
            data.interval
        )

        data.obs = pd.DataFrame(obs, index=samples_with_genotype.index)
        data.exp = pd.DataFrame(exp, index=samples_with_genotype.index)
        data.fdr = pd.DataFrame(fdr, index=samples_with_genotype.index)
        data.w = pd.DataFrame(w, index=samples_with_genotype.index)
        
        if calc_posteriors:
            post = self.calc_posteriors(
                obs, exp, fdr, w, data.disp_models
            )
            post = pd.DataFrame(post, index=samples_with_genotype.index)
            data.z = 1 - np.exp(-post)
        return data
    
    def load_fp_data(self, tabix_files, interval: GenomicInterval):
        """
        Analogous to PosteriorStats._load_data but implemented with TabixExtractor
        """
        n = len(tabix_files)
        m = len(interval)

        obs = np.zeros((n, m), dtype=np.float64)
        exp = np.zeros((n, m), dtype=np.float64)
        fdr = np.ones((n, m), dtype=np.float64)
        w = np.zeros((n, m), dtype=np.float64)


        for tabix_path in tabix_files:
            with TabixExtractor(
                tabix_path,
                columns=[
                    'chrom', 'start', 'end', 
                    'exp', 'obs', 'neglog_pval', 
                    'neglog_winpval', 'fdr'
                ]
            ) as extractor:
                tabix_data = extractor[interval]
                interval_index = tabix_data['start'].values - interval.start
                obs[:, interval_index] = tabix_data['obs']
                exp[:, interval_index] = tabix_data['exp']
                fdr[:, interval_index] = tabix_data['fdr']
                w[:, interval_index] = 1

        return obs, exp, fdr, w

    def calc_posteriors(self, obs, exp, fdr, w, disp_models):
        prior = posterior.compute_prior_weighted(fdr, w, cutoff=0.05) 
        scale = posterior.compute_delta_prior(
            obs, exp, fdr, 
            np.ones((obs.shape[0], 2)),
            cutoff=0.05
        )
        ll_on = posterior.log_likelihood(
            obs, exp, disp_models,
            delta=scale, w=3
        ) 
        ll_off = posterior.log_likelihood(
            obs, exp, disp_models,
            w=3
        )
        post = -posterior.posterior(prior, ll_on, ll_off)
        post[post <= 0] = 0.0
        return post
    

class DifferentialFootprintLoader(PlotDataLoader):

    def _load(self, data: DataBundle):

        # Store number of samples
        L_a = data.a_allele_count

        obs = np.ascontiguousarray(data.obs)
        exp = np.ascontiguousarray(data.exp)
        
        log2_obs_over_exp = np.log2((obs + 1) / (exp + 1)) # sample_id x positions
        
        #filter outliers per groups
        outliers = np.apply_along_axis(
            self.find_outliers,
            axis=0,
            arr=log2_obs_over_exp
        ) # sample_id x positions
        log2_obs_over_exp[outliers] = np.nan

        # Step 1: Fit prior
        variance = np.nanvar(log2_obs_over_exp, axis=0) # per position
        variance = variance[np.isfinite(variance) & (variance > 0)]

        
        obj = lambda p: -invchi2.log_likelihood(variance, p[0], p[1])
        nu_0, sig2_0 = scipy.optimize.fmin(obj, x0=[1.0, 1.0], disp=False)
        
        step_args = (-6, 6, 100)
        # Negative binomial probabilities
        nb = differential.compute_logpmf_values(
            data.disp_models,
            obs,
            exp,
            *step_args
        )
        
        # Calculate prior
        pr_a = differential.compute_log_prior_t(
            log2_obs_over_exp[:L_a, :],
            nu_0, sig2_0, 
            *step_args
        )
        pr_b = differential.compute_log_prior_t(
            log2_obs_over_exp[L_a:, :],
            nu_0, sig2_0, 
            *step_args
        )
        pr_ab = differential.compute_log_prior_t(
            log2_obs_over_exp, 
            nu_0, sig2_0, 
            *step_args
        )

        x = np.linspace(*step_args)
        mua = x[np.argmax(pr_a, axis=0)]
        mub = x[np.argmax(pr_b, axis=0)]

        data.lfc = mub - mua
        
        # psuedo-integration over 'depletion' scores
        pa = pr_a[:, :, np.newaxis] + nb[:, :, :L_a]
        pb = pr_b[:, :, np.newaxis] + nb[:, :, L_a:]
        pab = pr_ab[:, :, np.newaxis] + nb[:, :, :]

        # likelihood
        La = np.sum(logsumexp(pa, axis=0), axis=1)
        Lb = np.sum(logsumexp(pb, axis=0), axis=1)
        Lab = np.sum(logsumexp(pab, axis=0), axis=1)
        llr = La + Lb - Lab

        lrt_pvalue = st.chi2.sf(2 * llr, df=3)

         # include 1 in stouffer's aggregation
        lrt_pvalue = np.clip(lrt_pvalue, None, 1.0 - 1e-6)
        data.neglog10_pval = -np.log10(stouffers_z(lrt_pvalue, 3))

        return data

    @staticmethod
    def median_abs_dev(x):
        med = np.median(x)
        return med, np.median(np.abs(x - med))

    def find_outliers(self, x):
        med, mad = self.median_abs_dev(x)
        if mad > 0:
            dev = 0.6745 * (x - med) / mad
            return np.abs(dev) > 3.5
        else:
            return np.zeros(len(x), dtype=bool)


class AnnotateWithDHSLoader(PlotDataLoader):
    # Not even used... move somewhere else?
    def _load(self, data: DataBundle, dhs_anndata: ad.AnnData, max_n=1):
        # Density: Get density based on DHSs
        dhs_idx = filter_df_to_interval(dhs_anndata.var).index
        # Currently choosing DHSs on random, could do closest(?)
        if max_n is not None:
            sampled_indices = np.random.choice(
                len(dhs_idx),
                size=min(len(dhs_idx), max_n),
                replace=False,
                seed=42
            )
            dhs_idx = dhs_idx[sampled_indices]
        
        data.dhs_density = pd.DataFrame(
            dhs_anndata[:, dhs_idx].layers["density"],
            index=dhs_anndata.obs.index,
            columns=dhs_idx
        )
        return data
