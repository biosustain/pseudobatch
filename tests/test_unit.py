import pytest

from pseudobatch.data_correction import pseudobatch_transform


def test_input_contain_nan(simulated_fedbatch):
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
