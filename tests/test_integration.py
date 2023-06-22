## This scripts contain integration tests of all fedbatch data correction function.
import logging

import numpy as np
import pytest
import statsmodels.api as sm
from patsy import dmatrices
import pandas as pd

from pseudobatch.data_correction import (
    pseudobatch_transform,
    pseudobatch_transform_pandas,
)
from pseudobatch.datasets import load_standard_fedbatch


def fit_ols_model(formula_like: str, data: pd.DataFrame) -> sm.regression.linear_model.RegressionResultsWrapper:
    y, X = dmatrices(formula_like, data)
    model = sm.OLS(endog=y, exog=X)
    res = model.fit()
    return res

###################### INTEGRATION TESTS ##################################
def test_pseudobatch_transform_full_data_set(simulated_fedbatch):
    # correct glucose data
    simulated_fedbatch["corrected_glucose"] = pseudobatch_transform(
        measured_concentration=simulated_fedbatch["c_Glucose"].to_numpy(),
        reactor_volume=simulated_fedbatch["v_volume_before_sample"].to_numpy(),
        accumulated_feed=simulated_fedbatch["v_feed_accum"].to_numpy(),
        concentration_in_feed=93.75,
        sample_volume=simulated_fedbatch["sample_volume"]
        .fillna(0)
        .to_numpy(),  # the sample volume column contains nan when at times where no sample was taken
    )

    # correct biomass data
    simulated_fedbatch["corrected_biomass"] = pseudobatch_transform(
        measured_concentration=simulated_fedbatch["c_Biomass"].to_numpy(),
        reactor_volume=simulated_fedbatch["v_volume_before_sample"].to_numpy(),
        accumulated_feed=simulated_fedbatch["v_feed_accum"].to_numpy(),
        concentration_in_feed=0,
        sample_volume=simulated_fedbatch["sample_volume"]
        .fillna(0)
        .to_numpy(),  # the sample volume column contains nan when at times where no sample was taken
    )

    ## Calculate growth rate
    model_dat = simulated_fedbatch
    y, X = dmatrices(
        formula_like="np.log(corrected_biomass) ~ timestamp", data=model_dat
    )
    model_corrected = sm.OLS(endog=y, exog=X)
    res_corrected = model_corrected.fit()
    mu_hat = res_corrected.params[1]

    ## Calculate glucose yield coefficient
    model_dat = simulated_fedbatch
    y, X = dmatrices(
        formula_like="corrected_glucose ~ corrected_biomass", data=model_dat
    )
    model_corrected = sm.OLS(endog=y, exog=X)
    res_corrected = model_corrected.fit()
    Yxs = res_corrected.params[1]

    assert mu_hat == pytest.approx(
        0.1, 0.005
    )  # mu_hat = 0.1004353052206084, the non exact result is due to slightly changing actual growth rate in simulation. This comes from the monod kinetics of the growth rate.
    assert Yxs == pytest.approx(-3.70, 0.01)


def test_pseudobatch_transform_measurements_only(
    simulated_fedbatch_measurements_only,
):
    # correct glucose data
    simulated_fedbatch_measurements_only[
        "corrected_glucose"
    ] = pseudobatch_transform(
        measured_concentration=simulated_fedbatch_measurements_only[
            "c_Glucose"
        ].to_numpy(),
        reactor_volume=simulated_fedbatch_measurements_only[
            "v_volume_before_sample"
        ].to_numpy(),
        accumulated_feed=simulated_fedbatch_measurements_only[
            "v_feed_accum"
        ].to_numpy(),
        concentration_in_feed=93.75,
        sample_volume=simulated_fedbatch_measurements_only["sample_volume"]
        .fillna(0)
        .to_numpy(),  # the sample volume column contains nan when at times where no sample was taken
    )

    # correct biomass data
    simulated_fedbatch_measurements_only[
        "corrected_biomass"
    ] = pseudobatch_transform(
        measured_concentration=simulated_fedbatch_measurements_only[
            "c_Biomass"
        ].to_numpy(),
        reactor_volume=simulated_fedbatch_measurements_only[
            "v_volume_before_sample"
        ].to_numpy(),
        accumulated_feed=simulated_fedbatch_measurements_only[
            "v_feed_accum"
        ].to_numpy(),
        concentration_in_feed=0,
        sample_volume=simulated_fedbatch_measurements_only["sample_volume"]
        .fillna(0)
        .to_numpy(),  # the sample volume column contains nan when at times where no sample was taken
    )

    ## Calculate growth rate
    model_dat = simulated_fedbatch_measurements_only
    y, X = dmatrices(
        formula_like="np.log(corrected_biomass) ~ timestamp", data=model_dat
    )
    model_corrected = sm.OLS(endog=y, exog=X)
    res_corrected = model_corrected.fit()
    mu_hat = res_corrected.params[1]

    ## Calculate glucose yield coefficient
    model_dat = simulated_fedbatch_measurements_only
    y, X = dmatrices(
        formula_like="corrected_glucose ~ corrected_biomass", data=model_dat
    )
    model_corrected = sm.OLS(endog=y, exog=X)
    res_corrected = model_corrected.fit()
    Yxs = res_corrected.params[1]

    assert mu_hat == pytest.approx(
        0.1, 0.006
    )  # mu_hat = 0.1004353052206084, the non exact result is due to slightly changing actual growth rate in simulation. This comes from the monod kinetics of the growth rate.
    assert Yxs == pytest.approx(-3.70, 0.01)


def test_pseudobatch_transform_multiple_feeds(simulated_multiple_feeds):
    """
    Test if the correct Yxs can be retrieved from the pseudo batch transformed data from a simulation
    that contain multiple feeds.
    """
    glucose_in_feed1 = simulated_multiple_feeds["c_Glucose_feed1"].iloc[0]
    glucose_in_feed2 = 0

    pseudo_df = pseudobatch_transform_pandas(
        simulated_multiple_feeds.fillna({"sample_volume": 0}),
        ["c_Biomass", "c_Glucose"],
        "v_Volume",
        ["v_Feed_accum1", "v_Feed_accum2"],
        [[0, 0], [glucose_in_feed1, glucose_in_feed2]],
        "sample_volume",
    )

    ## Calculate the consumed glucose
    pseudo_df["c_Glucose_consumed_pseudo"] = (
        pseudo_df["c_Glucose_pseudo"].iloc[0] - pseudo_df["c_Glucose_pseudo"]
    )

    ## Calculate glucose yield coefficient
    model_dat = pseudo_df
    y, X = dmatrices(
        formula_like="c_Glucose_consumed_pseudo ~ c_Biomass_pseudo",
        data=model_dat,
    )
    model_corrected = sm.OLS(endog=y, exog=X)
    res_corrected = model_corrected.fit()
    Yxs = res_corrected.params[1]

    assert Yxs == pytest.approx(simulated_multiple_feeds["Yxs"].iloc[0])

def test_accepts_nan():
    """Test if the pseudobatch_transform function handles nan values correctly, 
    i.e. calculates the mu_hat when the nan are present in the data. This mimicks
    the situation where not all species where measured in all samples."""

    fedbatch_df = (load_standard_fedbatch()
        .query("sample_volume > 0")
        .reset_index(drop=True)
        # change every 2nd biomass measurement to nan
        .assign(c_Biomass = lambda df: pd.Series([
            x if idx % 2 != 0 else np.nan for idx, x in df['c_Biomass'].items()
        ]))
    )
    logging.debug(fedbatch_df.filter(['timestamp', 'c_Biomass', 'v_Volume']))
    fedbatch_df['c_Biomass_pseudo'] = pseudobatch_transform(
        fedbatch_df['c_Biomass'].to_numpy(),
        fedbatch_df['v_Volume'].to_numpy(),
        fedbatch_df['v_Feed_accum'].to_numpy(),
        0,
        fedbatch_df['sample_volume'].to_numpy(),
    )
    
    fedbatch_df['c_Glucose_pseudo'] = pseudobatch_transform(
        fedbatch_df['c_Glucose'].to_numpy(),
        fedbatch_df['v_Volume'].to_numpy(),
        fedbatch_df['v_Feed_accum'].to_numpy(),
        fedbatch_df.s_f.iloc[0],
        fedbatch_df['sample_volume'].to_numpy(),
    )

    growth_rate_model = fit_ols_model(
        "np.log(c_Biomass_pseudo) ~ timestamp",
        fedbatch_df,
    )

    substrate_yield_model = fit_ols_model(
        "c_Glucose_pseudo ~ c_Biomass_pseudo",
        fedbatch_df,
    )

    # fetching the true values used in the simulation
    Yxs_true = fedbatch_df.Yxs.iloc[0]
    mu_true = fedbatch_df.mu0.iloc[0]

    assert growth_rate_model.params[1] == pytest.approx(mu_true, 1e-6) 
    assert np.abs(substrate_yield_model.params[1]) == pytest.approx(Yxs_true, 1e-6)
     
