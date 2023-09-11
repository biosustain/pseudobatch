import pytest
import pandas as pd
import numpy as np

from pseudobatch import (
    pseudobatch_transform,
    pseudobatch_transform_pandas,
    metabolised_amount,
    hypothetical_concentration,
)
from pseudobatch.datasets import (
    load_standard_fedbatch,
    load_product_inhibited_fedbatch,
    load_cho_cell_like_fedbatch,
    load_real_world_yeast_fedbatch,
    load_volatile_compounds_fedbatch,
)
import logging
logging.basicConfig(level=logging.DEBUG)


def test_load_standard_fedbatch_unique_timestamps():
    df = load_standard_fedbatch()
    assert df.empty is False
    assert df["timestamp"].duplicated().sum() == 0


def test_load_product_inhibited_fedbatch_unique_timestamps():
    df = load_product_inhibited_fedbatch()
    assert df.empty is False
    assert df["timestamp"].duplicated().sum() == 0


def test_load_cho_cell_like_fedbatch_unique_timestamps():
    df = load_cho_cell_like_fedbatch()
    assert df.empty is False
    assert df["timestamp"].duplicated().sum() == 0


def test_load_real_world_yeast_fedbatch():
    """Test that the dataset is loaded correctly and that the shape is correct."""
    df = load_real_world_yeast_fedbatch()
    logging.debug(df.shape)
    assert df.empty is False
    assert df.shape == (11400, 12), "The dataset has changed. Update the test."


def test_load_volatile_compounds_fedbatch_unique_timestamps():
    df = load_volatile_compounds_fedbatch()
    assert df.empty is False
    assert df["timestamp"].duplicated().sum() == 0


def test_pseudobatch_transform_pandas_preserves_index():
    """Test that the index of the input dataframe is preserved in the output dataframe."""
    df = load_standard_fedbatch()

    # change index to something else
    df.index = np.arange(1000, 1000 + len(df))
    transformed_df = pseudobatch_transform_pandas(
        df=df,
        measured_concentration_colnames="c_Glucose",
        reactor_volume_colname="v_Volume",
        accumulated_feed_colname="v_Feed_accum",
        concentration_in_feed=df.s_f.iloc[0],
        sample_volume_colname="sample_volume",
    )
    assert df.index.equals(transformed_df.index)


def test_pseudobatch_transform_pandas_validation_missing_concentration_in_feed():
    """Test that the validation fails when the number of measured_concentration_colnames is
    not equal to the number of length of concentration in feed."""
    df = load_standard_fedbatch()

    # missing concentration in feed data
    with pytest.raises(ValueError) as _:
        pseudobatch_transform_pandas(
            df=df,
            measured_concentration_colnames=["c_Glucose", 'c_Biomass'],
            reactor_volume_colname="v_Volume",
            accumulated_feed_colname="v_Feed_accum",
            concentration_in_feed=df.s_f.iloc[0],
            sample_volume_colname="sample_volume",
        )


def test_pseudobatch_transform_pandas_validation_multiple_feeds():
    """Tests that validation fails if the concentration_in_feed is incorectly formatted,
    when multiple feeds are used.
    
    In this test the inner lists in concentration_in_feed iterates over the measured_concentration_colnames
    this is WRONG. The inner lists should iterate over the feeds."""

    df = load_cho_cell_like_fedbatch()

    with pytest.raises(ValueError) as _:
        pseudobatch_transform_pandas(
            df=df,
            measured_concentration_colnames=["c_Glucose", 'c_Biomass', "c_Glutamine"],
            reactor_volume_colname="v_Volume",
            accumulated_feed_colname=["v_Feed_accum", "v_Feed_accum_2"],
            concentration_in_feed=[[df.c_Glucose_feed1.iloc[0] , 0, df.c_Glutamine_feed1], [0, 0, df.c_Glutamine_feed2]],
            sample_volume_colname="sample_volume",
        )
