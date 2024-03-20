import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numba import njit
import warnings


@njit
def sample_multinomial(n_array, p_matrix, num_samples, seed=0, use_seed=True):
    """
    Sample for each n and row of ps with parallel execution.
    
    Args:
        n_array (np.ndarray): Array of counts to sample.
        p_matrix (np.ndarray): Matrix of probabilities to sample from.
        num_samples (int): Number of samples to generate.
        
    Returns:
        np.ndarray: Matrix of sample counts.  
    """
    num_categories = p_matrix.shape[1]
    results = np.zeros((num_samples, num_categories), dtype=np.int32)

    for i in range(num_samples):
        for j in range(len(n_array)):
            n = n_array[j]
            p = p_matrix[j]
            counts = np.zeros(len(p), dtype=np.int32)

            # Inlined categorical sampling
            for _ in range(n):
                if use_seed:
                    np.random.seed(seed)
                    seed += 10000
                r = np.random.random()
                cumulative = 0.0
                for k, transition_p in enumerate(p):
                    cumulative += transition_p
                    if r < cumulative:
                        counts[k] += 1
                        break

            results[i] += counts
    return results

@njit
def stratified_sampling_from_counts(counts_to_sample, all_counts, with_replacement=False, total_indices=None, seed=0, use_seed=True):
    """
    Perform stratified sampling with or without replacement, using a provided total index array or default sequential indices,
    with the sampling process repeated n_samples times. Each sampling iteration uses a corresponding row from counts_to_sample.

    Parameters:
    - counts_to_sample: np.array, shape (n_samples, n_bins) - The counts to sample for each sampling iteration, across bins.
    - all_counts: np.array, shape (n_bins,) - The total counts available in each bin.
    - with_replacement: bool - Whether sampling should be done with replacement.
    - total_indices: np.array, optional - A single array of indices for all records. Defaults to sequential indices based on all_counts.
    - seed: int - Seed for the random number generator.
    - use_seed: bool - Whether to use a seed for random number generation.

    Returns:
    - samples: np.array, shape (sum(counts_to_sample[0, :]), n_samples) - Each column corresponds to one iteration of sampling,
      containing indices of the sampled records. Each row within a column represents a sampled index from one of the bins.
    """
    if total_indices is None:
        total_indices = np.arange(all_counts.sum())

    n_samples = counts_to_sample.shape[0]
    total_samples = np.sum(counts_to_sample[0, :])
    samples = np.zeros((total_samples, n_samples), dtype=np.int32)
    bin_starts = np.cumsum(np.append(np.array([0], dtype=np.int32), all_counts))[:-1]

    for sample_iter in range(n_samples):
        position = 0
        for bin_index in range(all_counts.shape[0]):
            count = counts_to_sample[sample_iter, bin_index]
            start_idx = bin_starts[bin_index]
            end_idx = start_idx + all_counts[bin_index]
            bin_indices = total_indices[start_idx:end_idx]
            if use_seed:
                np.random.seed(seed)
                seed += 10000
            sampled_indices = np.random.choice(bin_indices, size=count, replace=with_replacement)
            samples[position:position+count, sample_iter] = sampled_indices
            position += count

    return samples


def get_probs(bin_counts, w):
    """
    Get transition probabilities based on the distance between bins.
    Probabilities are computed using a Gaussian kernel with width w.
    
    Args:
        bin_distances (np.ndarray): Matrix of distances between bins.
        w (float): Width of the Gaussian kernel to use for the distance-based probabilities.
    
    Returns:
        np.ndarray: Matrix of probabilities.
    """
    bin_distances = squareform(pdist(np.array(bin_counts.index.tolist()), 'euclidean'))
    transformed = np.exp(-bin_distances / (2 * w**2))
    return transformed / transformed.sum(axis=1)


def perturb_bin_counts(bin_counts, w, num_samples, seed=0, use_seed=True):
    """
    Perturb bin counts by sampling from a multinomial distribution with probabilities
    based on the distance between bins.

    Args:
        bin_counts (pd.Series): Bin counts to perturb.
        w (float): Width of the Gaussian kernel to use for the distance-based probabilities.
        num_samples (int): Number of perturbed samples to generate.

    Returns:
        np.ndarray: Matrix of perturbed bin counts.

    Shape:
        - Input: (n,)
        - Output: (num_samples, n)
    """
    if w > 0:
        probs = get_probs(bin_counts, w)
        return sample_multinomial(bin_counts.values, probs, num_samples, seed, use_seed)
    else:
        return np.tile(bin_counts.values, (num_samples, 1))
    

@njit
def indices_to_indicators(original_index, sampled_indices):
    """
    Convert sampled indices to a count-based indicator matrix with Numba's njit decorator,
    where values represent the count of times an index is sampled.
    
    Parameters:
    - original_index: np.array, shape (len(original_index),) - Array of all elements in order.
    - sampled_indices: np.array, shape (N, num_samples) - Indices of sampled elements.
    
    Returns:
    - sampled_indicators: np.array, shape (len(original_index), num_samples) - Count-based indicator matrix.
    """
    num_samples = sampled_indices.shape[1]
    len_original = len(original_index)
    
    sampled_indicators = np.zeros((len_original, num_samples), dtype=np.int32)
    
    for sample_num in range(num_samples):
        for idx in range(sampled_indices.shape[0]):
            sampled_indicators[sampled_indices[idx, sample_num], sample_num] += 1
    
    return sampled_indicators


def stratified_sampling(sampling_data, ref_data, matching_fields, num_samples=100, w=0, starting_seed=None, input_sorted=False,
                return_indicators=False, replace=False, n_to_sample=None, ignore_missing_bins=False):
    """
    Sample from sampling_data to match the distribution of ref_data.
    
    Args:
        sampling_data (pd.DataFrame): Data to sample from.
        ref_data (pd.DataFrame): Data to match the distribution of.
        matching_fields (list): Fields to match the distribution on.
        num_samples (int): Number of samples to generate.
        w (float): Width of the Gaussian kernel to use for the distance-based probabilities.
        starting_seed (int): Random seed. For each sample seed is incremented by 10000.
        input_sorted (bool): Whether the input data is already sorted by matching_fields.
        return_indicators (bool): Whether to return a matrix of sample indicators.
        replace (bool): Whether to sample with replacement.
        n_to_sample (int): Number of records to sample. If None, defaults to the total number of records in ref_data.
    Returns:
        np.ndarray: Matrix of sample indicators.
    """
    original_index = sampling_data.index.to_numpy()
    sorted_variants = sampling_data.reset_index(drop=True)
    if not input_sorted:
        sorted_variants = sorted_variants.sort_values(matching_fields)
        original_index = original_index[sorted_variants.index.to_numpy()]

    reference_bin_counts = ref_data[matching_fields].value_counts().sort_index()
    sampling_bin_counts = sorted_variants[matching_fields].value_counts().sort_index()

    total_reference = reference_bin_counts.sum()

    # make sure both reference and sampling bins are the same, fill with zeros
    all_bin_indices = reference_bin_counts.index.union(sampling_bin_counts.index).sort_values()
    reference_bin_counts = reference_bin_counts.reindex(all_bin_indices, fill_value=0)
    sampling_bin_counts = sampling_bin_counts.reindex(all_bin_indices, fill_value=0)

    not_in_sampling = (reference_bin_counts > 0) & (sampling_bin_counts == 0)
    if not_in_sampling.any():
        message = f"""Reference data contains {len(not_in_sampling[not_in_sampling].index.tolist())} bins not present in sampling data;
                      in total {reference_bin_counts[not_in_sampling].sum()} records ({reference_bin_counts[not_in_sampling].sum() / total_reference * 100:.1f}%)."""
        if ignore_missing_bins:
            warnings.warn(message + " These bins will be ignored when sampling data.")
        else:
            raise ValueError(message + " To ignore these bins, set ignore_missing_bins=True.")
        
    reference_bin_counts[not_in_sampling] = 0
    
    # workaround to pass to numba function
    use_seed = starting_seed is not None
    starting_seed = starting_seed if use_seed else 0

    if n_to_sample is None:
        n_to_sample = total_reference

    if replace:
        reference_bin_counts = np.round(reference_bin_counts * (n_to_sample / total_reference)).astype(int)
    else:
        assert n_to_sample % total_reference == 0, f"Number of records to sample ({n_to_sample}) must be a multiple of the total reference records ({total_reference}) when sampling without replacement."
        reference_bin_counts = reference_bin_counts * (n_to_sample // total_reference)
        assert (sampling_bin_counts >= reference_bin_counts).all(), "Sampling data must have at least as many records as reference data in each bin when sampling without replacement."

    bin_counts_to_sample = perturb_bin_counts(reference_bin_counts, w=w, num_samples=num_samples, seed=starting_seed, use_seed=use_seed)
    if not replace:
        assert (bin_counts_to_sample >= reference_bin_counts.values).all(), "After applying the perturbation, the number of records to sample from each bin must be less than or equal to the number of records in the sampling data."

    sampled = stratified_sampling_from_counts(
        bin_counts_to_sample, sampling_bin_counts.values, 
        with_replacement=replace,
        total_indices=original_index,
        seed=starting_seed,
        use_seed=use_seed
    )

    if return_indicators:
        return indices_to_indicators(original_index, sampled)
    else:
        return sampled
