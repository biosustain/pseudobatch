from enum import Enum
from importlib.resources import files
from typing import Dict, List, Optional, Tuple, Union

import arviz as az
import numpy as np

from cmdstanpy import CmdStanModel
from numpy.typing import NDArray
from pandas.core.generic import is_number
from pydantic import BaseModel, Field, root_validator

from pseudobatch import stan
from pseudobatch.util import (
    get_lognormal_params_from_quantiles,
    get_normal_params_from_quantiles,
)

KNOWN_QUANTITIES = {
    "sigma_v": 0.05,
    "sigma_c": [0.05] * 3,
    "sigma_f": 0.05,
    "sigma_s": 0.05,
    "sigma_cfeed": 0.05,
}
STAN_FILE = files(stan).joinpath("error_propagation.stan")


class Distribution0d(str, Enum):
    normal = "normal"
    lognormal = "lognormal"


class Prior0d(BaseModel):
    distribution: Distribution0d
    loc: Optional[float] = None
    pct1: Optional[float] = None
    pct99: Optional[float] = None
    scale: Optional[float] = Field(gt=0)

    @root_validator
    def check_locscale_or_pcts(cls, values):
        """Either loc and scale should be specified or pct1 and pct99."""
        missing_locscale = [f for f in ["loc", "scale"] if values[f] is None]
        missing_pcts = [f for f in ["pct1", "pct99"] if values[f] is None]
        assert len(missing_locscale) != 1, f"Missing {missing_locscale[0]}"
        assert len(missing_pcts) != 1, f"Missing {missing_pcts[0]}"
        assert not (
            (len(missing_locscale) == 0) and (len(missing_pcts) == 0)
        ), "Prior input is all None."
        return values

    def __init__(self, **data):
        super().__init__(**data)
        if self.loc is None:
            if self.distribution is Distribution0d.normal:
                assert self.pct1 is not None
                assert self.pct99 is not None
                self.loc, self.scale = get_normal_params_from_quantiles(
                    p1=0.01, x1=self.pct1, p2=0.99, x2=self.pct99
                )
            elif self.distribution is Distribution0d.lognormal:
                assert self.pct1 is not None
                assert self.pct99 is not None
                self.loc, self.scale = get_lognormal_params_from_quantiles(
                    p1=0.01, x1=self.pct1, p2=0.99, x2=self.pct99
                )


class Prior0dNormal(Prior0d):
    distribution = Distribution0d.normal


class Prior0dLogNormal(Prior0d):
    distribution = Distribution0d.lognormal


class PriorInput(BaseModel):
    prior_apump: Prior0dNormal
    prior_as: Prior0dNormal
    prior_v0: Prior0dLogNormal
    prior_f_nonzero: Prior0dLogNormal
    prior_cfeed_nonzero: Optional[Prior0dLogNormal] = None


def run_error_propagation(
    y_concentration: NDArray,
    y_reactor_volume: NDArray,
    y_feed_in_interval: NDArray,
    y_sample_volume: NDArray,
    sd_concentration: List[float],
    sd_reactor_volume: float,
    sd_feed_in_interval: float,
    sd_sample_volume: float,
    prior_input: dict,
    y_concentration_in_feed: Optional[float] = None,
    sd_concentration_in_feed: Optional[float] = None,
    species_names: Optional[List[Union[str,int]]] = None,
) -> az.InferenceData:
    """Run the error propagation analysis, returning and InferenceData object.

    Parameters
    ----------

    y_concentration : The measured concentrations. Note that this should be
    a 2d array.

    y_reactor_volume : Array of reactor volume measurements.

    y_feed_in_interval : Array of accumulated feed measurements.

    y_sample_volume : Array of sample volume measurements.

    sd_concentration : List of concentration measurement errors.

    sd_reactor_volume : Reactor volume measurement error.

    sd_feed_in_interval : Feed in interval measurement_error.

    sd_sample_volume : Sample volume measurement error.

    prior_input : Dictionary that can be used to load a PriorInput object.

    y_concentration_in_feed : Optional single numerical measurement of feed
    concentration.

    sd_concentration_in_feed : Optional error for feed
    concentration measurements. 

    species_names: Optional List of species names. Must match the number of
    species with measured concentration.

    Returns
    -------
    az.InferenceData

    """
    N, S = y_concentration.shape
    pi = PriorInput.parse_obj(prior_input)
    if y_concentration_in_feed is not None:
        assert sd_concentration_in_feed is not None, (
            "sd_concentration_in_feed can only be None "
            "if y_concentration_in_feed is also None."
        )
        assert pi.prior_cfeed_nonzero is not None, (
            "prior_cfeed_nonzero can only be None "
            "if y_concentration_in_feed is also None."
        )
    else:
        y_concentration_in_feed = 0.
        sd_concentration_in_feed = 1.
    if pi.prior_cfeed_nonzero is None:
        prior_cfeed_nonzero = Prior0dLogNormal(loc=0, scale=1)
    else:
        prior_cfeed_nonzero = pi.prior_cfeed_nonzero
    prior_m = [Prior0dLogNormal(pct1=1e-9, pct99=1e9) for _ in range(S)]
    data = {
        "N": N,
        "S": S,
        "y_v": y_reactor_volume,
        "y_c": y_concentration,
        "y_s": y_sample_volume,
        "y_f": y_feed_in_interval,
        "y_cfeed": y_concentration_in_feed,
        "sigma_v": sd_reactor_volume,
        "sigma_c": sd_concentration,
        "sigma_s": sd_sample_volume,
        "sigma_f": sd_feed_in_interval,
        "sigma_cfeed": sd_concentration_in_feed,
        "prior_apump": [pi.prior_apump.loc, pi.prior_apump.scale],
        "prior_as": [pi.prior_as.loc, pi.prior_as.scale],
        "prior_v0": [pi.prior_v0.loc, pi.prior_v0.scale],
        "prior_m": [[p.loc for p in prior_m], [p.scale for p in prior_m]],
        "prior_f_nonzero": [pi.prior_f_nonzero.loc, pi.prior_f_nonzero.scale],
        "prior_cfeed_nonzero": [
            prior_cfeed_nonzero.loc,
            prior_cfeed_nonzero.scale,
        ],
    }
    if species_names is None:
        species_names = list(range(S))
    coords = {"sample": list(range(N)), "species": species_names}
    dims = {
        "m": ["sample", "species"],
        "as": ["sample"],
        "s": ["sample"],
        "f": ["sample"],
        "c": ["sample", "species"],
        "v": ["sample"],
        "pseudobatch_c": ["sample", "species"]
    }
    data_prior = {**data, **{"likelihood": 0}}
    data_posterior = {**data, **{"likelihood": 1}}
    model = CmdStanModel(stan_file=STAN_FILE)
    mcmc_prior = model.sample(data=data_prior, show_progress=False)
    mcmc_posterior = model.sample(data=data_posterior, show_progress=False)
    return az.from_cmdstanpy(mcmc_posterior, prior=mcmc_prior, coords=coords, dims=dims)
