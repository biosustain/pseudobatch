"""Run a script that does error propagation."""
import argparse
import os

import arviz as az
import pandas as pd
from cmdstanpy import CmdStanModel  # type: ignore
from scipy.special import logit

from matplotlib import pyplot as plt

from src.util import (
    get_lognormal_parameters_from_quantiles,
    get_normal_parameters_from_quantiles,
)

HERE = (
    os.path.dirname(os.path.abspath(__file__))
    if "__file__" in locals()
    else "."
)
FEDBATCH_FILE = os.path.join(HERE, "..", "tests", "test_data", "fed-batch3.csv")
IDATA_FILE = os.path.join(HERE, "idata.json")
# value for glucose below
# C_FEED_MEASURED = 93.75
C_FEED_MEASURED = 0
STAN_FILE = os.path.join(HERE, "..", "src", "stan", "model.stan")
DEFAULT_SAMPLE_KWARGS = {
    "chains": 4,
    "iter_warmup": 2000,
    "iter_sampling": 1000,
    "adapt_delta": 0.9999,
    "metric": "dense",
}
PRIORS = {
    # lognormal locations and scales
    "sigma_v": [-3, 1],  # quantiles (0.025, 0.5, 0.975) = (0.007, 0.05, 0.35)
    "sigma_f": [-3, 1],  # quantiles (0.025, 0.5, 0.975) = (0.007, 0.05, 0.35)
    "sigma_c": [-3, 1],  # quantiles (0.025, 0.5, 0.975) = (0.007, 0.05, 0.35)
    "sigma_s": [-3, 1],  # quantiles (0.025, 0.5, 0.975) = (0.007, 0.05, 0.35)
    "sigma_cfeed": [
        -3,
        1,
    ],  # quantiles (0.025, 0.5, 0.975) = (0.007, 0.05, 0.35)
    "pump_bias": [0, 0.02],  # quantiles (0.025, 0.5, 0.975) = (0.76, 1, 1.32)
    "s_frac_logit": list(
        get_normal_parameters_from_quantiles(
            logit(0.05), 0.01, logit(0.4), 0.99
        )
    ),
    "v0": list(get_lognormal_parameters_from_quantiles(750, 0.01, 850, 0.99)),
    "m": list(get_lognormal_parameters_from_quantiles(400, 0.01, 16000, 0.99)),
    "f_nonzero": list(
        get_lognormal_parameters_from_quantiles(10, 0.01, 1000, 0.99)
    ),
    "cfeed_nonzero": list(
        get_lognormal_parameters_from_quantiles(0.01, 0.01, 0.1, 0.99)
    ),
}


def main():
    """Run main function"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_kwargs", nargs="+")
    argv = parser.parse_args()
    if argv.sample_kwargs is not None:
        custom_sample_kwargs = {
            a.split("=")[0]: a.split("=")[1] for a in argv.sample_kwargs
        }
    else:
        custom_sample_kwargs = {}
    df = (
        pd.read_csv(FEDBATCH_FILE, index_col=0)
        .dropna(subset=["sample_volume"])
        .drop_duplicates(subset=["timestamp"], keep="first")
        .assign(
            v_feed_interval=lambda df: (
                df["v_feed_accum"]
                .subtract(df["v_feed_accum"].shift(1))
                .fillna(0)
            )
        )
    )
    data = {
        "N": len(df),
        "y_v": df["v_volume"].values.tolist(),
        "y_f": df["v_feed_interval"].values.tolist(),
        "y_c": df["m_Biomass"].divide(df["v_volume"]).values.tolist(),
        "y_s": df["sample_volume"].values.tolist(),
        "y_cfeed": C_FEED_MEASURED,
        "t": df["timestamp"].values.tolist(),
        "prior_sigma_v": PRIORS["sigma_v"],
        "prior_sigma_f": PRIORS["sigma_f"],
        "prior_sigma_c": PRIORS["sigma_c"],
        "prior_sigma_s": PRIORS["sigma_s"],
        "prior_sigma_cfeed": PRIORS["sigma_cfeed"],
        "prior_pump_bias": PRIORS["pump_bias"],
        "prior_s_frac_logit": PRIORS["s_frac_logit"],
        "prior_v0": PRIORS["v0"],
        "prior_m": PRIORS["m"],
        "prior_f_nonzero": PRIORS["f_nonzero"],
        "prior_cfeed_nonzero": PRIORS["cfeed_nonzero"],
        "likelihood": 1,
    }
    nonzero_time_steps = [t for t, f in zip(data["t"], data["y_f"]) if f > 0]
    print(nonzero_time_steps)
    model = CmdStanModel(stan_file=STAN_FILE)
    sample_kwargs = {**DEFAULT_SAMPLE_KWARGS, **custom_sample_kwargs}
    mcmc = model.sample(data=data, **sample_kwargs)
    idata = az.from_cmdstanpy(
        mcmc,
        coords={
            "time": df["timestamp"].values.tolist(),
            # "nonzero": nonzero_time_steps,
        },
        dims={
            "v": ["time"],
            "f": ["time"],
            "c": ["time"],
            "s": ["time"],
            "s_frac_logit": ["time"],
            # "f_nonzero": ["nonzero"],
            "m": ["time"],
            "pseudobatch_c": ["time"],
        },
        observed_data={"v": data["y_v"], "f": data["y_f"], "c": data["y_c"]},
    )
    idata.to_json(IDATA_FILE)
    df_out = df.copy().assign(
        c_Glucose=lambda df: df["m_Glucose"] / df["v_volume"],
        pseudobatch_mean=idata.posterior["pseudobatch_c"]
        .mean(dim=["chain", "draw"])
        .values,
    )
    f = az.plot_pair(idata, var_names=["v", "f"])
    plt.gcf().savefig("fig.png")
    print(df_out.round(2))


if __name__ == "__main__":
    main()
