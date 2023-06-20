from enum import Enum
from importlib.resources import files
from typing import Dict, List, Optional, Union

import arviz as az
from cmdstanpy import CmdStanModel
from numpy.typing import NDArray
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
    prior_alpha_pump: Prior0dNormal
    prior_alpha_s: Prior0dNormal
    prior_v0: Prior0dLogNormal
    prior_m: List[Prior0dLogNormal]
    prior_f_nonzero: Prior0dLogNormal
    prior_cfeed_nonzero: Prior0dLogNormal


def run_error_propagation(
    measured_concentration: NDArray,
    reactor_volume: NDArray,
    accumulated_feed: NDArray,
    concentration_in_feed: Union[NDArray, float],
    sample_volume: NDArray,
    prior_input: dict,
    known_quantities: Dict = KNOWN_QUANTITIES,
) -> az.InferenceData:
    """Run the error propagation analysis, returning and InferenceData object.

    Parameters
    ----------
    measured_concentration : NDArray

    reactor_volume : NDArray

    accumulated_feed : NDArray

    concentration_in_feed : Union[NDArray, float]

    sample_volume : NDArray

    prior_input : PriorInput

    known_quantities : Dict

    Returns
    -------
    az.InferenceData

    """
    N, S = measured_concentration.shape
    pi = PriorInput.parse_obj(prior_input)
    if not isinstance(concentration_in_feed, float):
        raise ValueError(
            "Error propagation currently only supports a single feed concentration:"
            " please make sure this input is a float not an NDArray."
        )
    data = {
        "N": N,
        "S": S,
        "y_v": reactor_volume,
        "y_c": measured_concentration,
        "y_s": sample_volume,
        "y_f": accumulated_feed,
        "y_cfeed": concentration_in_feed,
        "sigma_v": known_quantities["sigma_v"],
        "sigma_c": known_quantities["sigma_c"],
        "sigma_s": known_quantities["sigma_s"],
        "sigma_f": known_quantities["sigma_f"],
        "sigma_cfeed": known_quantities["sigma_cfeed"],
        "prior_alpha_pump": [
            pi.prior_alpha_pump.loc,
            pi.prior_alpha_pump.scale,
        ],
        "prior_alpha_s": [pi.prior_alpha_s.loc, pi.prior_alpha_s.scale],
        "prior_v0": [pi.prior_v0.loc, pi.prior_v0.scale],
        "prior_m": [[p.loc for p in pi.prior_m], [p.scale for p in pi.prior_m]],
        "prior_f_nonzero": [pi.prior_f_nonzero.loc, pi.prior_f_nonzero.scale],
        "prior_cfeed_nonzero": [
            pi.prior_cfeed_nonzero.loc,
            pi.prior_cfeed_nonzero.scale,
        ],
        "likelihood": 1,
    }
    model = CmdStanModel(stan_file=STAN_FILE)
    mcmc = model.sample(data=data)
    return az.from_cmdstanpy(mcmc)
