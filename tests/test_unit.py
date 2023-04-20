import pytest
import pandas as pd
import numpy as np

from pseudobatch import (
    pseudobatch_transform,
    pseudobatch_transform_pandas,
)
from pseudobatch.datasets import (
    load_standard_fedbatch,
    load_product_inhibited_fedbatch,
    load_cho_cell_like_fedbatch,
)


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

def test_pseudobatch_transform_pandas_preserves_index():
    """Test that the index of the input dataframe is preserved in the output dataframe."""
    df = load_standard_fedbatch()

    # change index to something else
    df.index = np.arange(1000, 1000 + len(df))
    transformed_df = pseudobatch_transform_pandas(
        df=df,
        measured_concentration_colnames=["c_Glucose"],
        reactor_volume_colname="v_Volume",
        accumulated_feed_colname="v_Feed_accum",
        concentration_in_feed=[df.s_f.iloc[0]],
        sample_volume_colname="sample_volume",
    )
    assert df.index.equals(transformed_df.index)