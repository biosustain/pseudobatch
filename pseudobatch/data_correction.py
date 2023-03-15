import pandas as pd
import numpy as np
from typing import Union, Iterable
from numpy.typing import NDArray

def mass_removal_correction(
    c_meas: NDArray,
    v_samples: NDArray,
    v_meas: int,
    v_init: int,
    v_feed: int,
    c_feed: int,
):
    """Correct for the removal of mass."""
    current_conc = c_meas.iloc[-1]

    return (
        current_conc * v_meas - v_feed * c_feed + np.sum(c_meas * v_samples)
    ) / v_init


def mass_removal_correction_pandas_df(
    df,
    c_meas_colname: str,
    v_samples_colname: str,
    v_meas_colname: str,
    v_feed_colname: str,
    c_feed: int,
) -> NDArray:
    """Apply the pseudo_batch_transform function to a pandas dataframe.

    Expanding windows on dataframes where tricky to make work, therefore this
    manual implementation.
    """
    corrected_conc = np.array([])
    first_idx = df.index[0]
    for idx, _ in df.iterrows():
        corrected_conc = np.append(
            corrected_conc,
            mass_removal_correction(
                c_meas=df.loc[:idx, c_meas_colname],
                v_samples=df.loc[:idx, v_samples_colname],
                v_meas=df.loc[idx, v_meas_colname],
                v_feed=df.loc[idx, v_feed_colname],
                v_init=df.loc[first_idx, v_meas_colname]
                + df.loc[first_idx, v_samples_colname],
                # In the Timothy's implementation, he adds the t=0 measurement.
                # I want to avoid that.
                c_feed=c_feed,
            ),
        )

    return corrected_conc


def dilution_factor_correction(c_meas, adf, v_feed, c_feed, v_meas):
    return c_meas * adf - ((v_feed * c_feed) / v_meas) * adf


def shift(xs, n):
    """Shifts the elements of a numpy array by n steps.

    Numpy version of pandas .shift method.

    Origin:
    https://stackoverflow.com/questions/30399534/shift-elements-in-a-numpy-array
    """
    if n >= 0:
        return np.concatenate((np.full(n, np.nan), xs[:-n]))
    else:
        return np.concatenate((xs[-n:], np.full(-n, np.nan)))


def accumulated_dilution_factor(
    after_sample_reactor_volume: NDArray, sample_volume: NDArray
) -> NDArray:
    """Calculates the accumulated dilution factor.

    after_sample_reactor_volume is the volume of the bioreactor. At the sampling
    time points, this has to be the volume AFTER the sample was taken.
    """

    dilution_factor = (sample_volume + after_sample_reactor_volume) / shift(
        after_sample_reactor_volume, 1
    )
    dilution_factor = np.nan_to_num(
        dilution_factor, copy=True, nan=1
    )  # The dilution factor of the first index is set to 1
    return np.cumprod(dilution_factor)


def pseudo_batch_transform(
    measured_concentration: NDArray,
    reactor_volume: NDArray,
    accumulated_feed: NDArray,
    concentration_in_feed: Union[NDArray, float],
    sample_volume: NDArray,
) -> NDArray:
    """Does the pseudo batch transformation for one species.

    The function returns a NDArray of the pseudo batch transformed
    measurements.

    Parameters
    ----------
    measured_concentration : NDArray

        a NDArray with the measured concentration of the species that should be
        transformed, e.g. biomass or compound

    reactor_volume : NDArray

        a NDArray of bioreator volume. The volume MUST be the volumes just
        BEFORE sampling.

    accumulated_feed : NDArray

        a NDArray of the accumulated volumen of feed added at this timepoint.

    concentration_in_feed : Union[NDArray, float]

        a NDArray OR a float of the concentration of the species in the feed,
        for biomass and products this is 0

    sample_volume : NDArray

        a NDArray of the sample volumes at given time points. The array should
    contain 0 at timepoints where no samples was taken

    """
    after_sample_reactor_volume = reactor_volume - sample_volume

    for i in [
        measured_concentration,
        reactor_volume,
        accumulated_feed,
        sample_volume,
    ]:
        if np.isnan(i).sum() > 0:
            msg = (
                "Nan was found in input data. Replace nan with an appropriate"
                " number - this is often 0."
            )
            raise ValueError(msg)
        adf = accumulated_dilution_factor(
            after_sample_reactor_volume, sample_volume
        )

    def fed_species_term(
        accum_feed: NDArray,
        conc_in_feed: Union[NDArray, float],
        reactor_vol: NDArray,
    ) -> NDArray[np.float64]:
        """Calculate the feed in interval.

        Prepand 0 to make the length of the array the same as the other arrays.

        """
        feed_in_interval = np.diff(accum_feed, prepend=0)
        return adf * feed_in_interval * conc_in_feed / reactor_vol

    if len(np.shape(accumulated_feed)) == 1:
        return measured_concentration * adf - np.cumsum(
            fed_species_term(
                accumulated_feed, concentration_in_feed, reactor_volume
            )
        )

    # Iterate over the columns of the accumulated feed to calculate
    # fed_species_term for each feed. Then calculate the row-wise sum of the
    # fed_species_term.

    fed_species_term_collection = np.zeros(np.shape(accumulated_feed)[0])
    for i in range(np.shape(accumulated_feed)[1]):
        fed_species_term_collection += fed_species_term(
            accumulated_feed[:, i], concentration_in_feed[i], reactor_volume
        )

    return measured_concentration * adf - np.cumsum(fed_species_term_collection)


def pseudo_batch_transform_multiple(
    measured_concentrations: NDArray,
    reactor_volume: NDArray,
    accumulated_feed: NDArray,
    concentration_in_feed: Union[NDArray, NDArray],
    sample_volume: NDArray,
) -> NDArray:
    """Does the pseudo batch transformation for several species.

    measured_concetration:

    a NDArray with the measured concentration of the species that should be
    transformed, e.g. biomass or compound. Each column should be a different
    species.

    reactor_volume:

    a NDArray of bioreator volume. The volume MUST be the volumes just BEFORE
    sampling.

    feed_in_interval:

    a NDArray of the feed volumen in the interval since last timepoint.

    concentration_in_feed:

    a NDArray with the concentrations of the species in the feed. The order has
    to match the order of the species given in measured_concentration.
    Alternatively an NDArray wich contain the concetration for each
    individual time step. Again the order of the columns has to match in order
    in measured_concentration.

    sample_volume:

    a NDArray of the sample volumes at given time points. The array should
    contain 0 at timepoints where no samples was taken
    """
    bad_shape_msg = (
        "concentration_in_feed needs to be 2D. Either a row vector or column"
        " vector was given. Try using .reshape(1,-1)"
    )
    assert len(concentration_in_feed.shape) == 2, bad_shape_msg
    out_list = list()
    for col in range(0, measured_concentrations.shape[1]):
        transformed = pseudo_batch_transform(
            measured_concentration=measured_concentrations[:, col],
            reactor_volume=reactor_volume,
            accumulated_feed=accumulated_feed,
            concentration_in_feed=concentration_in_feed[:, col],
            sample_volume=sample_volume,
        )
        out_list.append(np.array([transformed]))
    return np.concatenate(out_list).transpose()


def pseudo_batch_transform_pandas(
    df: pd.DataFrame,
    measured_concentration_colnames: Union[str, Iterable[str]],
    reactor_volume_colname: str,
    accumulated_feed_colname: Union[str, Iterable[str]],
    concentration_in_feed: Union[Iterable[float], Iterable[NDArray[np.float64]]],
    sample_volume_colname: str,
    pseudo_col_postfix: str = "_pseudo",
) -> pd.DataFrame:
    """Apply pseudo batch transformation for several species from a dataframe.

    The function outputs a dataframe containing the pseudo batch transformed
    data a species for each column.

    Parameters
    ----------
    df : pd.DataFrame

        a pandas dataframe containing the data to be transformed

    measured_concentration_colnames : Union[str, Iterable[str]]

        a string or list of strings with the column names of the measured
        concentration of the species that should be transformed, e.g. biomass or
        glucose

    reactor_volume_colname : str

        a string with the column name of the reactor volume

    accumulated_feed_colname : Union[str, Iterable[str]]

        a string or list of strings with the column names of the accumulated feed

    concentration_in_feed : Union[Iterable[float], Iterable[Iterable[float]]]

        a list of floats or list of lists of floats with the concentrations of
        the species in the feed. The order has to match the order of the species
        given in measured_concentration_colnames. If multiple feed streams are
        given, the concentration_in_feed should be a list of lists. E.g.
        [[conc_Glc_feed1, conc_Glc_feed2], [conc_Biomass_feed1,
        conc_Biomass_feed2]]. Thus, the outer list iterates over the feed
        streams and the inner list iterates over the species. The order has to
        match the order of the species, and accumulated feeds given in
        measured_concentration_colnames and accumulated_feed_colname.

    sample_volume_colname : str

        a string with the column name of the sample volume

    pseudo_col_postfix : str, optional

        a string with the postfix to be added to the column names of the pseudo
        batch transformed data, by default "_pseudo"

    Returns
    -------

    pd.DataFrame

        a pandas dataframe containing the pseudo batch transformed data a
        species for each column.

    """
    out = pd.DataFrame()
    for species, conc in zip(
        measured_concentration_colnames, concentration_in_feed
    ):
        out[species + pseudo_col_postfix] = pseudo_batch_transform(
            measured_concentration=df[species].to_numpy(),
            reactor_volume=df[reactor_volume_colname].to_numpy(),
            accumulated_feed=df[accumulated_feed_colname].to_numpy(),
            concentration_in_feed=conc,
            sample_volume=df[sample_volume_colname].to_numpy(),
        )
    return out
