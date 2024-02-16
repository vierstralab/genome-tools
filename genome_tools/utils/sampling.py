import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from numba import jit, prange


@jit(nopython=True, parallel=False)
def sample_multinomial(n_array, p_matrix, num_samples):
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

    for i in prange(num_samples):
        for j in range(len(n_array)):
            n = n_array[j]
            p = p_matrix[j]
            counts = np.zeros(len(p), dtype=np.int32)

            # Inlined categorical sampling
            for _ in range(n):
                r = np.random.random()
                cumulative = 0.0
                for k, p_val in enumerate(p):
                    cumulative += p_val
                    if r < cumulative:
                        counts[k] += 1
                        break

            results[i] += counts
    return results

@jit(nopython=True, parallel=False)
def get_sample_indicators(counts_to_sample, all_counts, seed=0):
    """
    Sample from counts_to_sample to match the distribution of all_counts.
    
    Args:
        counts_to_sample (np.ndarray): Matrix of counts to sample from.
        all_counts (np.ndarray): Array of total counts for each column.
        seed (int): Random seed.
    
    Returns:
        np.ndarray: Matrix of sample indicators.
    """
    total_rows = all_counts.sum()
    num_samples, num_columns = counts_to_sample.shape
    res = np.zeros((total_rows, num_samples), dtype=np.bool_)

    cums = 0
    for i in range(num_columns):
        all_c = all_counts[i]
        for j in range(num_samples):
            np.random.seed(seed)
            c_to_sample = counts_to_sample[j, i]
            binary_matrix = np.arange(all_c) < c_to_sample
            shuffled_indices = np.random.permutation(all_c)
            res[cums: cums + all_c, j] = binary_matrix[shuffled_indices]
            seed += 10000
        cums += all_c

    return res



def get_probs(bin_distances, w):
    """
    Get probabilities based on the distance between bins.
    
    Args:
        bin_distances (np.ndarray): Matrix of distances between bins.
        w (float): Width of the Gaussian kernel to use for the distance-based probabilities.
    
    Returns:
        np.ndarray: Matrix of probabilities.
    """
    transformed = np.exp(-bin_distances / (2*w**2))
    return transformed / transformed.sum(axis=1)

def fast_sample(sampling_data, ref_data, matching_fields, num_samples=100, w=0, input_sorted=False):
    """
    Sample from sampling_data to match the distribution of ref_data.
    
    Args:
        sampling_data (pd.DataFrame): Data to sample from.
        ref_data (pd.DataFrame): Data to match the distribution of.
        matching_fields (list): Fields to match the distribution on.
        num_samples (int): Number of samples to generate.
        w (float): Width of the Gaussian kernel to use for the distance-based probabilities.
        input_sorted (bool): Whether the input data is already sorted by matching_fields.
    
    Returns:
        np.ndarray: Matrix of sample indicators.
    """
    if input_sorted:
        sorted_variants = sampling_data.reset_index(drop=True)
        reordering_indices = np.arange(sampling_data.shape[0])
    else:
        sorted_variants = sampling_data.reset_index(drop=True).sort_values(matching_fields)
        reordering_indices = np.argsort(sorted_variants.index.to_numpy()) 

    reference_bin_counts = ref_data[matching_fields].value_counts().sort_index()
    all_bin_counts = sorted_variants[matching_fields].value_counts().sort_index()
    bin_counts_to_sample = perturb_bin_counts(reference_bin_counts, w=w, num_samples=num_samples)
    
    sample_indicators = get_sample_indicators(bin_counts_to_sample, all_bin_counts.values).astype(bool)
    return sample_indicators[reordering_indices, :]

def perturb_bin_counts(bin_counts, w=0.01, num_samples=1000):
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
    bin_distances = squareform(pdist(np.array(bin_counts.index.tolist()), 'euclidean'))
    if w > 0:
        probs = get_probs(bin_distances, w)
        return sample_multinomial(bin_counts.values, probs, num_samples)
    else:
        return np.tile(bin_counts.values, (num_samples, 1))




def get_cumulative_dist(values):
    z  = np.sort(values)
    x = np.arange(z.shape[0])/z.shape[0]
    return z, x
