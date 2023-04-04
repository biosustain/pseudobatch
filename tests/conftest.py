# This file defines the fixtures for the pytest workflow
import pytest
import sys
import os
import pandas as pd
import numpy as np

sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir))
)
from pseudobatch.data_correction import shift
from pseudobatch.datasets import load_cho_cell_like_fedbatch


###################### SETUP SIMULATED DATA FIXTURES ################################
@pytest.fixture(scope="session")
def simulated_fedbatch():
    """Loads simulated fed-batch dataset and takes only a few points to mimick
    samples of true fermentation."""
    fedbatch_file = os.path.join("tests", "test_data", "fed-batch3.csv")
    fedbatch_df = pd.read_csv(fedbatch_file, index_col=0).drop_duplicates(
        subset="timestamp", keep="last"
    )  # ODE solver save both data before and after sampling event
    assert isinstance(fedbatch_df, pd.DataFrame)
    # Calculating the concentrations
    fedbatch_df["c_Biomass"] = (
        fedbatch_df["m_Biomass"] / fedbatch_df["v_volume"]
    )
    fedbatch_df["c_Glucose"] = (
        fedbatch_df["m_Glucose"] / fedbatch_df["v_volume"]
    )

    fedbatch_df["v_volume_before_sample"] = fedbatch_df[
        "v_volume"
    ] + fedbatch_df["sample_volume"].fillna(0)

    return fedbatch_df


@pytest.fixture(scope="session")
def simulated_fedbatch_measurements_only():
    """Loads simulated fed-batch dataset and takes only a few points to mimick
    samples of true fermentation."""
    fedbatch_file = os.path.join(
        "tests", "test_data", "fed-batch3_measurements_only.csv"
    )
    fedbatch_df = pd.read_csv(fedbatch_file, index_col=0)

    # Calculating the concentrations
    fedbatch_df["c_Biomass"] = (
        fedbatch_df["m_Biomass"] / fedbatch_df["v_volume"]
    )
    fedbatch_df["c_Glucose"] = (
        fedbatch_df["m_Glucose"] / fedbatch_df["v_volume"]
    )

    fedbatch_df["v_volume_before_sample"] = fedbatch_df[
        "v_volume"
    ] + fedbatch_df["sample_volume"].fillna(0)

    return fedbatch_df


@pytest.fixture(scope="session")
def simulated_multiple_feeds():
    """ """
    return load_cho_cell_like_fedbatch()
