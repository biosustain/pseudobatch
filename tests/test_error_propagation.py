from importlib.resources import files

import numpy as np
import pandas as pd
from scipy.special import logit

from pseudobatch import run_error_propagation
from pseudobatch.datasets import data

KNOWN_QUANTITIES = {
    "sigma_v": 0.05,
    "sigma_c": [0.05] * 3,
    "sigma_f": 0.05,
    "sigma_s": 0.05,
    "sigma_cfeed": 0.05,
}

PRIORS_GOOD = {
    "prior_apump": {"pct1": np.log(1 - 0.1), "pct99": np.log(1 + 0.1)},
    "prior_as": {"pct1": logit(0.05), "pct99": logit(0.4)},
    "prior_v0": {"pct1": 1000, "pct99": 1030},  # Here is the bad bit!
    "prior_m": [
        {"pct1": 200, "pct99": 200000},
        {"pct1": 2, "pct99": 2000},
        {"pct1": 200, "pct99": 200000},
    ],
    "prior_f_nonzero": {"pct1": 10, "pct99": 1000},
    "prior_cfeed_nonzero": {"pct1": 0.01, "pct99": 0.1},
}


def test_error_propagation():
    data_path = files(data).joinpath("standard_fed-batch_process.csv")
    species = ["Product", "Glucose", "Biomass"]
    samples = (
        pd.read_csv(str(data_path), index_col=0)
        .dropna(subset=["sample_volume"])
        .drop_duplicates(subset=["timestamp"], keep="first")
        .assign(
            v_feed_interval=lambda df: (
                df["v_Feed_accum"] - df["v_Feed_accum"].shift(1)
            ).fillna(0),
        )
        .reset_index()
    )
    idata = run_error_propagation(
        measured_concentration=samples[["c_" + s for s in species]].values,
        reactor_volume=samples["v_Volume"].values,
        accumulated_feed=samples["v_Feed_accum"].values,
        concentration_in_feed=0.0,
        sample_volume=samples["sample_volume"].values,
        prior_input=PRIORS_GOOD,
        known_quantities=KNOWN_QUANTITIES,
    )
