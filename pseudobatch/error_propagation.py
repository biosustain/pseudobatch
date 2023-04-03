from enum import Enum
from typing import Dict, List, Optional, Union

import arviz as az
from cmdstanpy import CmdStanModel
from numpy.typing import NDArray
from pydantic import BaseModel, Field, root_validator

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


class Distribution0d(str, Enum):
    normal = "normal"
    lognormal = "lognormal"


class Prior0d(BaseModel):
    distribution: Distribution0d
    loc: Optional[float]
    pct1: Optional[float]
    pct99: Optional[float]
    scale: Optional[float] = Field(gt=0)

    @root_validator
    def check_locscale_or_pcts(cls, values):
        """Either loc and scale should be specified or pct1 and pct99."""
        missing_locscale = [
            f for f in ["loc", "scale"] if f not in cls.__fields_set__
        ]
        missing_pcts = [
            f for f in ["pct1", "pct99"] if f not in cls.__fields_set__
        ]
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
    def __init__(self, loc, pct1, pct99, scale):
        return super().__init__(
            distribution=Distribution0d.normal,
            loc=loc,
            pct1=pct1,
            pct99=pct99,
            scale=scale,
        )


class Prior0dLogNormal(Prior0d):
    def __init__(self, loc, pct1, pct99, scale):
        return super().__init__(
            distribution=Distribution0d.lognormal,
            loc=loc,
            pct1=pct1,
            pct99=pct99,
            scale=scale,
        )


class PriorInput(BaseModel):
    prior_alpha_pump: Prior0d
    prior_alpha_s: Prior0d
    prior_v0: Prior0d
    prior_m: List[Prior0d]
    prior_f_nonzero: Prior0d
    prior_cfeed_nonzero: Prior0d


class ErrorPropagationConfig(BaseModel):
    prior_input: PriorInput
    known_quantities: Dict = KNOWN_QUANTITIES


def run_error_propagation(
    measured_concentration: NDArray,
    reactor_volume: NDArray,
    accumulated_feed: NDArray,
    concentration_in_feed: Union[NDArray, float],
    sample_volume: NDArray,
    epc: ErrorPropagationConfig,
) -> az.InferenceData:
    """Run the error propagation analysis, returning and InferenceData object.

    Parameters
    ----------
    measured_concentration : NDArray

    reactor_volume : NDArray

    accumulated_feed : NDArray

    concentration_in_feed : Union[NDArray, float]

    sample_volume : NDArray

    epc : ErrorPropagationConfig

    Returns
    -------
    az.InferenceData

    """
    model = CmdStanModel(stan_file="stan/error_propagation.stan")
    N, S = measured_concentration.shape
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
        "sigma_v": KNOWN_QUANTITIES["sigma_v"],
        "sigma_c": KNOWN_QUANTITIES["sigma_c"],
        "sigma_s": KNOWN_QUANTITIES["sigma_s"],
        "sigma_f": KNOWN_QUANTITIES["sigma_f"],
        "sigma_cfeed": KNOWN_QUANTITIES["sigma_cfeed"],
        "prior_alpha_pump": [
            epc.prior_input.prior_alpha_pump.loc,
            epc.prior_input.prior_alpha_pump.scale,
        ],
        "prior_alpha_s": [
            epc.prior_input.prior_alpha_s.loc,
            epc.prior_input.prior_alpha_s.scale,
        ],
        "prior_v0": [
            epc.prior_input.prior_v0.loc,
            epc.prior_input.prior_v0.scale,
        ],
        "prior_m": [
            [p.loc for p in epc.prior_input.prior_m],
            [p.scale for p in epc.prior_input.prior_m],
        ],
        "prior_f_nonzero": [
            epc.prior_input.prior_f_nonzero.loc,
            epc.prior_input.prior_f_nonzero.scale,
        ],
        "prior_cfeed_nonzero": [
            epc.prior_input.prior_cfeed_nonzero.loc,
            epc.prior_input.prior_cfeed_nonzero.scale,
        ],
        "likelihood": 1,
    }
    mcmc = model.sample(data=data)
    return az.from_cmdstanpy(mcmc)
