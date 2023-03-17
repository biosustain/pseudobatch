"""This modules provides functions to load the simulated datasets that are used in the tests and examples."""
import pandas as pd
import pathlib


def _load_simulated_dataset(filename: str) -> pd.DataFrame:
    """Load and prepare the simulated dataset. At the sampling time the 
    simulation contains a value both before and after the sample was taken.
    When the the data is loaded only the first value is kept. This is the value 
    before the sample was taken."""

    data_path = pathlib.Path(__file__).parent / "data" / filename
    fedbatch_df = (
        pd.read_csv(data_path)
        .drop_duplicates(subset="timestamp", keep="first")
        .reset_index(
            drop=True
        )  # Reset the index to avoid gaps in index numbers from removed rows.
        .fillna(
            {"sample_volume": 0}
        )  # Fill the sample volume column with 0 when no sample was taken.
    )
    return fedbatch_df


def load_standard_fedbatch():
    """Load the standard fed-batch process dataset. This dataset
    mimicks a substrate limited exponential fed-batch process utilizing
    a glucose feed. During the fed-batch process, samples are withdrawn.
    The parameters values use for the simulation is stored in the
    dataframe.
    """
    return _load_simulated_dataset("standard_fed-batch_process.csv")


def load_product_inhibeted_fedbatch():
    """Load the product inhibition fed-batch process dataset. This dataset
    mimicks a substrate limited exponential fed-batch process utilizing
    a glucose feed. In this simulation the cell growth is inhibited by
    increasing product concentration. During the fed-batch process,
    samples are withdrawn. The parameters values use for the simulation
    is stored in the dataframe.
    """
    return _load_simulated_dataset("product_inhibition.csv")


def load_cho_cell_like_fedbatch():
    """Load the CHO cell like fed-batch process dataset. This dataset
    mimicks fed-batch process carried out in a AMBR15 cultivation system.
    The main characteristic is that this simulation utilises two substrates
    both of which are fed. They are fed through two different feed media.
    During the fed-batch process, samples are withdrawn. The parameters
    values use for the simulation is stored in the dataframe.
    """
    return _load_simulated_dataset("multiple_impulse_feed_process.csv")
