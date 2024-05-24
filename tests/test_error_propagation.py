from importlib.resources import files

import numpy as np
import pandas as pd
from scipy.special import logit

from pseudobatch.error_propagation import run_error_propagation
from pseudobatch.datasets import load_standard_fedbatch


EXAMPLE_PRIOR_INPUT = {
    "prior_apump": {"pct1": np.log(1 - 0.1), "pct99": np.log(1 + 0.1)},
    "prior_as": {"pct1": logit(0.05), "pct99": logit(0.4)},
    "prior_v0": {"pct1": 1000, "pct99": 1030},  # Here is the bad bit!
    "prior_cfeed": [
        {"loc": 0, "scale": 1},
        {"pct1": 98, "pct99": 102},
        {"loc": 0, "scale": 1}
    ],
}


def test_error_propagation():
    
    species = ["Product", "Glucose", "Biomass"]
    samples = (
        load_standard_fedbatch(sampling_points_only=True)
        .assign(
            v_feed_interval=lambda df: (
                df["v_Feed_accum"] - df["v_Feed_accum"].shift(1)
            ).fillna(0),
        )
        .reset_index()
    )
    _ = run_error_propagation(
        y_concentration=samples[["c_" + s for s in species]].values,
        y_reactor_volume=samples["v_Volume"].values,
        y_feed_in_interval=samples["v_feed_interval"].values,
        y_concentration_in_feed=np.array([0., 100., 0.]),
        y_sample_volume=samples["sample_volume"].values,
        prior_input=EXAMPLE_PRIOR_INPUT,
        sd_reactor_volume=0.05,
        sd_concentration=[0.05] * 3,
        sd_feed_in_interval=0.05,
        sd_sample_volume=0.05,
        sd_concentration_in_feed=0.05,
        seed=123,
    )
