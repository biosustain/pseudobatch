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
    metabolised_amount,
    hypothetical_concentration,
)
from pseudobatch.datasets import (
    load_standard_fedbatch,
    load_product_inhibited_fedbatch,
    load_cho_cell_like_fedbatch,
    load_volatile_compounds_fedbatch,
    load_evaporation_fedbatch,
)


def fit_ols_model(formula_like: str, data: pd.DataFrame) -> sm.regression.linear_model.RegressionResultsWrapper:
    y, X = dmatrices(formula_like, data)
    model = sm.OLS(endog=y, exog=X)
    res = model.fit()
    return res

###################### INTEGRATION TESTS ##################################
def test_pseudobatch_transform_full_data_set():
    """Tests that the growth rate and substrate yields are correctly estimated 
    pseudo batch transformed data. This test utilises all simulated data points"""
    # correct glucose data
    fedbatch_df = load_standard_fedbatch() 
    fedbatch_df["corrected_glucose"] = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Glucose"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=fedbatch_df.s_f.iloc[0], # the glucose concentration in the feed is stored in the dataframe
        sample_volume=fedbatch_df["sample_volume"].to_numpy(),
    )

    # correct biomass data
    fedbatch_df["corrected_biomass"] = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Biomass"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=0,
        sample_volume=fedbatch_df["sample_volume"].to_numpy(),  
    )

    ## Calculate growth rate
    res_corrected = fit_ols_model(
        formula_like="np.log(corrected_biomass) ~ timestamp",
        data=fedbatch_df
    )
    mu_hat = res_corrected.params[1]

    ## Calculate glucose yield coefficient
    model_dat = fedbatch_df
    res_corrected = fit_ols_model(
        formula_like="corrected_glucose ~ corrected_biomass", data=model_dat
    )
    Yxs = np.abs(res_corrected.params[1]) # yields have to be positive

    # True values
    mu_true = fedbatch_df.mu_true.iloc[0]
    Yxs_true = fedbatch_df.Yxs.iloc[0]

    assert mu_hat == pytest.approx(
        mu_true, 1e-6
    ) 
    assert Yxs == pytest.approx(Yxs_true, 1e-6)


def test_pseudobatch_transform_sampling_points_only():
    fedbatch_df = load_standard_fedbatch(sampling_points_only=True)
    # correct glucose data
    fedbatch_df["corrected_glucose"] = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Glucose"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=fedbatch_df.s_f.iloc[0], # the glucose concentration in the feed is stored in the dataframe
        sample_volume=fedbatch_df["sample_volume"].to_numpy(), 
    )

    # correct biomass data
    fedbatch_df["corrected_biomass"] = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Biomass"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=0,
        sample_volume=fedbatch_df["sample_volume"].to_numpy(),
    )

    ## Calculate growth rate
    res_corrected = fit_ols_model(
        formula_like="np.log(corrected_biomass) ~ timestamp",
        data=fedbatch_df
    )
    mu_hat = res_corrected.params[1]

    ## Calculate glucose yield coefficient
    model_dat = fedbatch_df
    res_corrected = fit_ols_model(
        formula_like="corrected_glucose ~ corrected_biomass", data=model_dat
    )
    Yxs = np.abs(res_corrected.params[1]) # yields have to be positive

    # True values
    mu_true = fedbatch_df.mu_true.iloc[0]
    Yxs_true = fedbatch_df.Yxs.iloc[0]

    # When we only use the sampling points the accuracy is lower
    assert mu_hat == pytest.approx(
        mu_true, 1e-3
    ) 
    assert Yxs == pytest.approx(Yxs_true, 1e-3)


def test_pseudobatch_transform_multiple_feeds():
    """
    Test if the correct Yxs can be retrieved from the pseudo batch transformed data from a simulation
    that contain multiple feeds.
    """
    fedbatch_df = load_cho_cell_like_fedbatch()
    glucose_in_feed1 = fedbatch_df["c_Glucose_feed1"].iloc[0]
    glucose_in_feed2 = 0

    fedbatch_df[['pseudo_Biomass', 'pseudo_Glucose']] = pseudobatch_transform_pandas(
        fedbatch_df,
        ["c_Biomass", "c_Glucose"],
        "v_Volume",
        ["v_Feed_accum1", "v_Feed_accum2"],
        [[0, 0], [glucose_in_feed1, glucose_in_feed2]],
        "sample_volume",
    )

    ## Calculate glucose yield coefficient
    res_corrected = fit_ols_model(
        formula_like="pseudo_Glucose ~ pseudo_Biomass",
        data=fedbatch_df,
    )
    Yxs = np.abs(res_corrected.params[1])

    assert Yxs == pytest.approx(fedbatch_df["Yxs"].iloc[0])

def test_accepts_nan():
    """Test if the pseudobatch_transform function handles nan values correctly, 
    i.e. calculates the mu_hat when the nan are present in the data. This mimicks
    the situation where not all species where measured in all samples."""

    fedbatch_df = (load_standard_fedbatch(sampling_points_only=True)
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

    assert growth_rate_model.params[1] == pytest.approx(mu_true, 1e-4) 
    assert np.abs(substrate_yield_model.params[1]) == pytest.approx(Yxs_true, 1e-4)


def test_calculation_of_gaseous_yield():
    '''Tests that the gaseous yield can be correctly estimated using the 
    hypothetical_concentration(), metabolised_amount and the pseudobatch_transform functions.'''
    fedbatch_df = load_volatile_compounds_fedbatch()

    # Preprocess the data
    fedbatch_df['m_O2_after_sample'] = fedbatch_df['m_O2'] - fedbatch_df['c_O2'] * fedbatch_df['sample_volume']
    fedbatch_df['m_O2_consumed'] = metabolised_amount(
        off_gas_amount=fedbatch_df['m_O2_gas'].to_numpy(),
        dissolved_amount_after_sampling=fedbatch_df['m_O2_after_sample'].to_numpy(),
        inlet_gas_amount=fedbatch_df['m_O2_in'].to_numpy(),
        sampled_amount=(fedbatch_df['c_O2'] * fedbatch_df['sample_volume']).cumsum().to_numpy(),
    )
    fedbatch_df['hypothetical_c_O2'] = hypothetical_concentration(
        metabolised_amount=fedbatch_df['m_O2_consumed'].to_numpy(),
        reactor_volume=fedbatch_df['v_Volume'].to_numpy(),
        sample_volume=fedbatch_df['sample_volume'].to_numpy(),
    )

    # Pseudo batch transform
    fedbatch_df[['pseudo_Biomass','pseudo_O2']] = pseudobatch_transform_pandas(
        df=fedbatch_df,
        measured_concentration_colnames=['c_Biomass','hypothetical_c_O2'],
        reactor_volume_colname='v_Volume',
        accumulated_feed_colname='v_Feed_accum',
        concentration_in_feed=[0,0],
        sample_volume_colname='sample_volume',
    )

    # Estimate CO2 yield
    res = fit_ols_model(
        "pseudo_O2 ~ pseudo_Biomass",
        fedbatch_df
    )

    Yxco2_true = fedbatch_df.Yxo2.iloc[0]
    assert np.abs(res.params[1]) == pytest.approx(Yxco2_true, 1e-6)


def test_pseudobatch_transformation_evaporation():
    """Tests that the growth rate and substrate yields are correctly estimated
    pseudo batch transformed data. This test utilises all simulated data points
    """
    # correct glucose data
    fedbatch_df = load_evaporation_fedbatch()
    fedbatch_df["corrected_glucose"] = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Glucose"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=fedbatch_df.s_f.iloc[
            0
        ],  # the glucose concentration in the feed is stored in the dataframe
        sample_volume=fedbatch_df["sample_volume"].to_numpy(),
    )

    # correct biomass data
    fedbatch_df["corrected_biomass"] = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Biomass"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=0,
        sample_volume=fedbatch_df["sample_volume"].to_numpy(),
    )

    ## Calculate growth rate
    res_corrected = fit_ols_model(
        formula_like="np.log(corrected_biomass) ~ timestamp", data=fedbatch_df
    )
    mu_hat = res_corrected.params[1]

    ## Calculate glucose yield coefficient
    model_dat = fedbatch_df
    res_corrected = fit_ols_model(
        formula_like="corrected_glucose ~ corrected_biomass", data=model_dat
    )
    Yxs = np.abs(res_corrected.params[1])  # yields have to be positive

    # True values
    mu_true = fedbatch_df.mu_true.iloc[0]
    Yxs_true = fedbatch_df.Yxs.iloc[0]

    assert mu_hat == pytest.approx(mu_true, 1e-4)
    assert Yxs == pytest.approx(Yxs_true, 1e-6)
