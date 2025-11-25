import numpy as np
from scipy.stats import binned_statistic_2d
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize


def nan_gaussian_filter(arr, sigma):
    """
    Apply a Gaussian filter to an array with NaNs in a way that ignores NaNs.
    
    Parameters:
      arr: 2D numpy array with NaNs for missing values.
      sigma: Standard deviation for Gaussian kernel.
      
    Returns:
      The smoothed array, with the same NaN positions as appropriate.
    """
    # Create a mask of valid values (1 for valid, 0 for NaN)
    mask = np.isfinite(arr).astype(float)
    
    # Replace NaNs in the array with 0
    arr_zeroed = np.nan_to_num(arr, nan=0.0)
    
    # Convolve both the array and the mask with the Gaussian filter
    arr_conv = gaussian_filter(arr_zeroed, sigma=sigma, mode='constant', cval=0.0)
    mask_conv = gaussian_filter(mask, sigma=sigma, mode='constant', cval=0.0)
    
    # Avoid division by zero: where the mask is zero, we want to retain NaN
    with np.errstate(divide='ignore', invalid='ignore'):
        result = arr_conv / mask_conv
    result[mask_conv == 0] = np.nan
    
    return result

def compute_weighted_stats(Z, counts):
    """
    Compute weighted mean and standard deviation on array with non-finite values.

    Parameters
    ----------
    Z : array_like
        Array of values. Can contain NaN or infinite entries, which are ignored.
    counts : array_like
        Array of non-negative weights (e.g., counts or frequencies) aligned with `Z`.
        Must be broadcastable to the shape of `Z`.

    Returns
    -------
    mu : float
        Weighted mean of finite entries in `Z`.
    sigma : float
        Weighted standard deviation of finite entries in `Z`.
    """
    mask = np.isfinite(Z)
    values = Z[mask]
    wts = counts[mask]
    mu = np.sum(values * wts) / np.sum(wts)
    sigma = np.sqrt(np.sum(wts * (values - mu)**2) / np.sum(wts))
    return mu, sigma


def plot_contours(x, y, z, mask_radius=1.5, sigma=50, levels=50, linewidths=0.2, lw_outer=0, bins=1500, cmap='solar_extra', vmin=None, vmax=None, z_score=False, ax=None, rasterized=False, **kwargs):
    
    """
    Plot smoothed contour maps of scalar values over a 2D embedding.

    This function bins `(x, y)` coordinates into a regular grid, computes the mean
    of associated values `z` per bin, applies Gaussian smoothing, masks regions
    without nearby data points, and plots contour and filled contour maps.

    Parameters
    ----------
    x, y : array_like, shape (n_samples,)
        Coordinates of the embedding (e.g., UMAP or t-SNE).
    z : array_like, shape (n_samples,)
        Scalar values associated with each coordinate.
    mask_radius : float, default=0.5
        Maximum distance from grid points to data points for inclusion.
        Grid cells farther than this radius from any data point are masked (set to NaN).
    sigma : float, default=8
        Standard deviation for the Gaussian kernel used in smoothing (in grid cells).
    levels : int or array_like, default=15
        Number or positions of contour levels passed to `matplotlib.contour`/`contourf`.
    linewidths : float, default=0.2
        Line width of inner contour lines.
    lw_outer : float or None, default=None
        Line width of the outer boundary contour. If None, defaults to `linewidths`.
    bins : int or [int, int], default=1500
        Number of bins along each axis for the 2D histogram grid.
    cmap : str or Colormap, default='solar_extra'
        Colormap for the filled contours.
    vmin, vmax : float, optional
        Value range for colormap normalization. If None, inferred from masked grid.
    z_score : bool, default=False
        If True, transform the masked smoothed grid values into Z-scores,
        weighted by counts in each bin.
    ax : matplotlib.axes.Axes, optional
        Axes on which to draw the plot. If None, a new figure and axes are created.
    **kwargs
        Additional keyword arguments passed to `matplotlib.contour` and `contourf`.

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axes containing the plot.
    sm : matplotlib.cm.ScalarMappable
        Scalar mappable object for creating a colorbar.

    Notes
    -----
    - The grid is expanded by 5% beyond the min/max of the input coordinates to avoid edge effects.
    - Gaussian smoothing is applied via `nan_gaussian_filter`, which preserves NaNs.
    - Regions without nearby data (determined by `mask_radius`) are excluded from
      the contour plot to prevent extrapolation artifacts.
    - If `z_score=True`, weighted mean and standard deviation are computed using
      `compute_weighted_stats`, and the smoothed values are standardized.
    """

    if lw_outer is None:
        lw_outer = linewidths

    # Assume embedding (N x 2 array) and values (N array) are defined
    # embedding = your UMAP coordinates, values = your scalar values

    # Compute original limits
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    
    # Compute 5% expansion
    dx = 0.05 * (x_max - x_min)
    dy = 0.05 * (y_max - y_min)
    
    # Define the expanded range
    range_val = ((x_min - dx, x_max + dx), (y_min - dy, y_max + dy))
    
    # Create the value grid using binned statistic (mean)
    Z, x_edges, y_edges, _ = binned_statistic_2d(
        x, y, z, statistic='mean', bins=bins,
        range=range_val
    )
    
    # Apply Gaussian smoothing
    Z_smooth = nan_gaussian_filter(Z, sigma=sigma)
    
    # Generate meshgrid for plotting (using bin centers)
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2
    X, Y = np.meshgrid(x_centers, y_centers, indexing='ij')    
    
    # Build a KDTree for scatter points
    tree = cKDTree(np.stack([x, y]).T)
    
    # Flatten the grid points into a list of (x,y) pairs
    grid_points = np.column_stack((X.ravel(), Y.ravel()))
    
    # For each grid point, find the distance to the nearest scatter point
    distances, _ = tree.query(grid_points, distance_upper_bound=mask_radius)
    
    # Reshape distances to grid shape and create mask:
    # A grid point is kept if the distance is less than mask_radius.
    distance_mask = distances.reshape(X.shape) <= mask_radius
    
    # Mask the smoothed grid: grid points without a nearby scatter point are set to NaN
    Z_masked = np.where(distance_mask, Z_smooth, np.nan)

    # Transform to Z-scores
    if z_score:
        counts, _, _, _ = binned_statistic_2d(x, y, z, statistic='count', bins=bins, range=range_val)
        counts_masked = np.where(distance_mask, counts, np.nan)
        mu, sig = compute_weighted_stats(Z_masked, counts_masked)
        Z_masked = (Z_masked - mu) / sig
    
    

    # Plot scatter points for context
    if ax is None:
        ax = plt.gca()
    # plt.scatter(X[distance_mask], Y[distance_mask], s=2, c='w', rasterized=True)

    vmin = vmin if vmin is not None else np.nanmin(Z_masked)
    vmax = vmax if vmax is not None else np.nanmax(Z_masked)
    
    # Plot only the contour lines over the scatter, with white background outside scatter
    if linewidths > 0:
        contour_lines = ax.contour(X, Y, Z_masked, levels=levels, linewidths=linewidths, colors='grey', **kwargs)
        if rasterized:
            for c in contour_lines.collections:
                c.set_rasterized(True)

    
    cf = ax.contourf(
        X, Y, Z_masked,
        levels=levels, cmap=cmap, 
        vmin=vmin, vmax=vmax,
        linestyles='solid',
        **kwargs
    )

    if rasterized:
        for c in cf.collections:
            c.set_rasterized(True)

    ax.axis('off')

    outline_contour = ax.contour(X, Y, distance_mask.astype(float), levels=[0.5], colors='grey', linestyles='solid', linewidths=lw_outer, **kwargs)
    if rasterized:
        for c in outline_contour.collections:
            c.set_rasterized(True)

    
    # Now, create a separate figure with a colorbar.
    # We'll use the min and max of the masked grid (ignoring NaNs)

    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # dummy array for the colorbar
    
    return ax, sm


