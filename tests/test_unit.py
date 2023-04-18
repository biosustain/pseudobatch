import pytest
import pandas as pd
import numpy as np
from pseudobatch import pseudobatch_transform, reverse_pseudobatch_transform
from pseudobatch.datasets import (
    load_standard_fedbatch,
    load_product_inhibited_fedbatch,
    load_cho_cell_like_fedbatch,
)
import logging
logging.basicConfig(level=logging.DEBUG)

def test_input_contain_nan(simulated_fedbatch: pd.DataFrame):
    # correct glucose data
    with pytest.raises(ValueError) as _:
        pseudobatch_transform(
            measured_concentration=simulated_fedbatch["c_Glucose"].to_numpy(),
            reactor_volume=simulated_fedbatch[
                "v_volume_before_sample"
            ].to_numpy(),
            accumulated_feed=simulated_fedbatch["v_feed_accum"].to_numpy(),
            concentration_in_feed=93.75,
            sample_volume=simulated_fedbatch[
                "sample_volume"
            ].to_numpy(),  # the sample volume column contains nan when at times where no sample was taken
        )


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

def test_reverse_pseudobatch_transform():
    """Test that the reverse pseudobatch transform is the inverse of the 
    pseudobatch transform."""

    fedbatch_df = load_standard_fedbatch().query("sample_volume > 0")
    c_Glucose_pseudo = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Glucose"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=fedbatch_df.s_f.iloc[0],
        sample_volume=fedbatch_df["sample_volume"].to_numpy(),
    )

    c_Glucose_reverse_pseudo = reverse_pseudobatch_transform(
        pseudo_concentration=c_Glucose_pseudo,
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=fedbatch_df.s_f.iloc[0],
        sample_volume=fedbatch_df["sample_volume"].to_numpy(),
    )

    assert np.allclose(fedbatch_df["c_Glucose"].to_numpy(), c_Glucose_reverse_pseudo, atol=1e-7)
