import numpy as np


def weighted_block_jackknife(theta, theta_js, m_js):
    """Compute the jackknife standard error and Z-score from a weighted block jackknife procedure.

    Args:
        theta (float): The observed genome-wide statistic computed over all blocks.
        theta_js (numpy.ndarray): Leave-one-out estimates, one per block with informative sites.
        m_js (numpy.ndarray): Block weights (number of segregating sites per block).

    Returns:
        result (tuple): A tuple of (g, sigma, Z) where g is the number of valid blocks,
            sigma is the jackknife standard error, and Z is the Z-score (theta / sigma).
            Returns (nan, nan, nan) when inference is not possible.

    Notes
    -----
    Implements the weighted delete-m jackknife from Patterson (2005), where the
    Z-score uses the original genome-wide estimate theta (not the bias-corrected
    theta_J).

    Returns (nan, nan, nan) when g <= 1 (insufficient blocks for a variance
    estimate), n == 0 (no observations), or any h_j <= 1 (a single block
    contains all observations, making the denominator h_j - 1 invalid).
    """
    # Identify blocks with finite, positive-weight leave-one-out estimates.
    is_valid = np.isfinite(theta_js) & np.isfinite(m_js) & (m_js > 0)

    # Retain only valid blocks for all subsequent computation.
    theta_js = theta_js[is_valid]
    m_js = m_js[is_valid]

    # Number of valid blocks.
    g = m_js.size

    # Total number of observations across all valid blocks.
    n = np.sum(m_js)

    # Inverse fractional block sizes: h_j = n / m_j.
    h_js = n / m_js

    # Guard: need at least two blocks, nonzero total, and no block containing all observations.
    if (g <= 1) or (n == 0) or np.any(h_js <= 1):
        return np.nan, np.nan, np.nan

    # Jackknifed estimate (Patterson 2005 Eq. 1): bias-corrected point estimate.
    theta_j_hat = (g * theta) - np.sum(((n - m_js) / n) * theta_js)

    # Pseudovalues (Patterson 2005 Eq. 3): one per block.
    tau_js = (h_js * theta) - ((h_js - 1) * theta_js)

    # Variance estimate (Patterson 2005 Eq. 4): average weighted squared deviation.
    sigma = np.sqrt((1 / g) * np.sum(((tau_js - theta_j_hat) ** 2) / (h_js - 1)))

    # Z-score: genome-wide estimate divided by its standard error.
    z_score = theta / sigma if (sigma != 0) and np.isfinite(sigma) else np.nan

    return g, sigma, z_score
