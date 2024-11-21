"""This modules provides functions to load the simulated datasets that are used in the tests and examples."""
import pandas as pd
import pathlib


def _prepare_simulated_dataset(data_path: str, sampling_points_only: bool = False) -> pd.DataFrame:
    """Load and prepare the simulated dataset. At the sampling time the 
    simulation contains a value both before and after the sample was taken.
    When the the data is loaded only the first value is kept. This is the value 
    before the sample was taken.
    
    Parameters
    ----------
    data_path : str
        Path to the csv file containing the simulated data.
    sampling_points_only : bool, optional
        If True, only the rows where a sample was taken is kept, by default False
    """

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

    if sampling_points_only:
        fedbatch_df = fedbatch_df[fedbatch_df["sample_volume"] > 0].reset_index(
            drop=True
        )

    return fedbatch_df


def load_standard_fedbatch(sampling_points_only: bool = False):
    """Load the standard fed-batch process dataset. This dataset
    mimicks a substrate limited exponential fed-batch process utilizing
    a glucose feed. During the fed-batch process, samples are withdrawn.
    The parameters values use for the simulation is stored in the
    dataframe.

    Parameters
    ----------
    sampling_points_only : bool, optional
        If True, only the rows where a sample was taken is kept, by default False
    """
    data_path = pathlib.Path(__file__).parent / "data" / "standard_fed-batch_process.csv"
    return _prepare_simulated_dataset(data_path, sampling_points_only=sampling_points_only)


def load_product_inhibited_fedbatch(sampling_points_only: bool = False):
    """Load the product inhibition fed-batch process dataset. This dataset
    mimicks a substrate limited exponential fed-batch process utilizing
    a glucose feed. In this simulation the cell growth is inhibited by
    increasing product concentration. During the fed-batch process,
    samples are withdrawn. The parameters values use for the simulation
    is stored in the dataframe.

    Parameters
    ----------
    sampling_points_only : bool, optional
        If True, only the rows where a sample was taken is kept, by default False
    """
    data_path = pathlib.Path(__file__).parent / "data" / "product_inhibition.csv"
    return _prepare_simulated_dataset(data_path, sampling_points_only=sampling_points_only)


def load_cho_cell_like_fedbatch(sampling_points_only: bool = False):
    """Load the CHO cell like fed-batch process dataset. This dataset
    mimicks fed-batch process carried out in a AMBR15 cultivation system.
    The main characteristic is that this simulation utilises two substrates
    both of which are fed. They are fed through two different feed media.
    During the fed-batch process, samples are withdrawn. The parameters
    values use for the simulation is stored in the dataframe.

    Parameters
    ----------
    sampling_points_only : bool, optional
        If True, only the rows where a sample was taken is kept, by default False
    """
    data_path = pathlib.Path(__file__).parent / "data" / "multiple_impulse_feed_process.csv"
    return _prepare_simulated_dataset(data_path, sampling_points_only=sampling_points_only)


def load_real_world_yeast_fedbatch():
    """Load the real world yeast fed-batch process dataset. This dataset
    is obtained from an experiment carried out in a biolector."""

    data_path = pathlib.Path(__file__).parent / "data" / "biolector_yeast_fedbatch.csv"
    return pd.read_csv(data_path)


def load_volatile_compounds_fedbatch(sampling_points_only: bool = False):
    """Load the volatile compounds fed-batch process dataset. This dataset
    mimicks a substrate limited exponential fed-batch process utilizing
    a glucose feed. During the fed-batch process, samples are withdrawn.
    The parameters values use for the simulation is stored in the
    dataframe.

    Parameters
    ----------
    sampling_points_only : bool, optional
        If True, only the rows where a sample was taken is kept, by default False
    """
    data_path = pathlib.Path(__file__).parent / "data" / "volatile_product.csv"
    return _prepare_simulated_dataset(data_path, sampling_points_only=sampling_points_only)


def load_evaporation_fedbatch(sampling_points_only: bool = False):
    """Load the evaporation fed-batch process dataset. This dataset
    mimicks a substrate limited exponential fed-batch process utilizing
    a glucose feed. During the fed-batch process, samples are withdrawn.
    The parameters values use for the simulation is stored in the
    dataframe.

    Parameters
    ----------
    sampling_points_only : bool, optional
        If True, only the rows where a sample was taken is kept, by default False
    """
    data_path = pathlib.Path(__file__).parent / "data" / "evaporation.csv"
    return _prepare_simulated_dataset(data_path, sampling_points_only=sampling_points_only)