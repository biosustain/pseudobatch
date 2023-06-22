from typing import Tuple

import numpy as np
from scipy.stats import norm


def get_normal_params_from_quantiles(
    x1: float, p1: float, x2: float, p2: float
) -> Tuple[float, float]:
    """Find parameters for a normal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ normal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.
    """
    ppf_low, ppf_high = (norm.ppf(x, scale=1) for x in [p1, p2])
    mu = (x1 * ppf_high - x2 * ppf_low) / (ppf_high - ppf_low)
    sigma = (x2 - x1) / (ppf_high - ppf_low)
    return mu, sigma


def get_lognormal_params_from_quantiles(
    x1: float, p1: float, x2: float, p2: float
) -> Tuple[float, float]:
    """Find parameters for a lognormal distribution from two quantiles.

    i.e. get mu and sigma such that if X ~ lognormal(mu, sigma), then pr(X <
    x1) = p1 and pr(X < x2) = p2.
    """
    return get_normal_params_from_quantiles(np.log(x1), p1, np.log(x2), p2)
