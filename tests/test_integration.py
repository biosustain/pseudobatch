## This scripts contain integration tests of all fedbatch data correction function.
import pytest

import numpy as np
from pseudobatch.data_correction import (
    pseudo_batch_transform,
    pseudo_batch_transform_pandas,
)
from patsy import dmatrices
import statsmodels.api as sm


###################### INTEGRATION TESTS ##################################
def test_pseudo_batch_transform_full_data_set(simulated_fedbatch):
    # correct glucose data
    simulated_fedbatch["corrected_glucose"] = pseudo_batch_transform(
        measured_concentration=simulated_fedbatch["c_Glucose"].to_numpy(),
        reactor_volume=simulated_fedbatch["v_volume_before_sample"].to_numpy(),
        accumulated_feed=simulated_fedbatch["v_feed_accum"].to_numpy(),
        concentration_in_feed=93.75,
        sample_volume=simulated_fedbatch["sample_volume"]
        .fillna(0)
        .to_numpy(),  # the sample volume column contains nan when at times where no sample was taken
    )

    # correct biomass data
    simulated_fedbatch["corrected_biomass"] = pseudo_batch_transform(
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


def test_pseudo_batch_transform_measurements_only(
    simulated_fedbatch_measurements_only,
):
    # correct glucose data
    simulated_fedbatch_measurements_only[
        "corrected_glucose"
    ] = pseudo_batch_transform(
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
    ] = pseudo_batch_transform(
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


def test_pseudo_batch_transform_multiple_feeds(simulated_multiple_feeds):
    """
    Test if the correct Yxs can be retrieved from the pseudo batch transformed data from a simulation
    that contain multiple feeds.
    """
    glucose_in_feed1 = simulated_multiple_feeds["c_Glucose_feed1"].iloc[0]
    glucose_in_feed2 = 0

    pseudo_df = pseudo_batch_transform_pandas(
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
