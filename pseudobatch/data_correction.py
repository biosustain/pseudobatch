from typing import Iterable, Union

import numpy as np
import pandas as pd
from numpy.typing import NDArray


def _shift(xs, n):
    """Shifts the elements of a numpy array by n steps.

    Numpy version of pandas .shift method.

    Origin:
    https://stackoverflow.com/questions/30399534/shift-elements-in-a-numpy-array

    Parameters
    ----------
    xs : NDArray
        The array to be shifted
    n : int
        The number of steps to shift. Positive values shift to the right,
        negative values shift to the left.

    Returns
    -------
    NDArray
        The shifted array. The first n elements are set to np.nan if n is
        positive, the last n elements are set to np.nan if n is negative.
    """
    if n >= 0:
        return np.concatenate((np.full(n, np.nan), xs[:-n]))
    else:
        return np.concatenate((xs[-n:], np.full(-n, np.nan)))


def accumulated_dilution_factor(
    after_sample_reactor_volume: NDArray, sample_volume: NDArray
) -> NDArray:
    """Calculates the accumulated dilution factor.

    Parameters
    ----------
    after_sample_reactor_volume : NDArray
        The volume of the bioreactor. At the sampling time points, this has to
        be the volume AFTER the sample was taken.
    sample_volume : NDArray
        The volume of the sample taken at each time point.

    Returns
    -------
    NDArray
        The accumulated dilution factor for each timepoint. The first value is
        always 1.
    """

    dilution_factor = (sample_volume + after_sample_reactor_volume) / _shift(
        after_sample_reactor_volume, 1
    )
    dilution_factor = np.nan_to_num(
        dilution_factor, copy=True, nan=1
    )  # The dilution factor of the first index is set to 1
    return np.cumprod(dilution_factor)


def pseudobatch_transform(
    measured_concentration: NDArray,
    reactor_volume: NDArray,
    accumulated_feed: NDArray,
    concentration_in_feed: Union[NDArray, float],
    sample_volume: NDArray,
) -> NDArray:
    """Pseudo batch transformation function for a single species. This function
    transforms the measured concentrations to the pseudo concentrations.

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

    Returns
    -------
    NDArray
        A NDArray with the pseudo concentrations of the species.

    """
    after_sample_reactor_volume = reactor_volume - sample_volume

    for i in [
        reactor_volume,
        accumulated_feed,
        sample_volume,
    ]:
        if np.isnan(i).sum() > 0:
            msg = (
                "Nan was found in either the reactor volume, accumulated feed or "
                "the sample volume. Replace nan with an appropriate value."
                "For example, if no sample was taken, the sample volume should be "
                "set to 0."
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

        Prepend the first value of the accumulated feed to the array because the
        is the feed added during the first interval.
        """
        feed_in_interval = np.diff(accum_feed, prepend=accum_feed[0])
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


def pseudobatch_transform_multiple(
    measured_concentrations: NDArray,
    reactor_volume: NDArray,
    accumulated_feed: NDArray,
    concentration_in_feed: Union[Iterable, NDArray],
    sample_volume: NDArray,
) -> NDArray:
    """Perform the pseudo batch transformation on multiple species at once. This
    function simply wraps the `pseudobatch_transform()` function.

    Parameters
    ----------
    measured_concentrations : NDArray
        a NDArray with the measured concentration of the species that should be
        transformed, e.g. biomass or compound. Each column should be a different
        species and the reach row a time point.
    reactor_volume : NDArray
        a NDArray of bioreator volume. The volume MUST be the volumes just BEFORE
        sampling.
    feed_in_interval : NDArray
        a NDArray of the feed volumen in the interval since last timepoint.
    concentration_in_feed : Iterable or NDArray
        the concentration of each species in the feed medium. This can be an iterable
        of the same length as the number of columns in the measured_concentrations array.
        Alternatively an NDArray which contain the concetration for each
        individual time step. Again the order of the columns has to match in order
        in measured_concentration.
    sample_volume : NDArray
        a NDArray of the sample volumes at given time points. The array should
        contain 0 at timepoints where no samples was taken

    Returns
    -------
    NDArray
        A NDArray of the same shape as the measured_concentration argument with the
        pseudo concentrations of the species. The columns correspond to the columns
        in measured_concentration.
    """
    bad_shape_msg = (
        "concentration_in_feed needs to be 2D. Either a row vector or column"
        " vector was given. Try using .reshape(1,-1)"
    )
    assert len(concentration_in_feed.shape) == 2, bad_shape_msg
    out_list = list()
    for col in range(0, measured_concentrations.shape[1]):
        transformed = pseudobatch_transform(
            measured_concentration=measured_concentrations[:, col],
            reactor_volume=reactor_volume,
            accumulated_feed=accumulated_feed,
            concentration_in_feed=concentration_in_feed[:, col],
            sample_volume=sample_volume,
        )
        out_list.append(np.array([transformed]))
    return np.concatenate(out_list).transpose()


def pseudobatch_transform_pandas(
    df: pd.DataFrame,
    measured_concentration_colnames: Union[str, Iterable[str]],
    reactor_volume_colname: str,
    accumulated_feed_colname: Union[str, Iterable[str]],
    concentration_in_feed: Union[
        Iterable[Union[float, int]], Iterable[Iterable[Union[float, int]]]
    ],
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
        conc_Biomass_feed2]]. Thus, the outer list iterates over the species
        and the inner list iterates over the feed streams. The order has to
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

    # convert input to list if only one column name is given
    if isinstance(measured_concentration_colnames, str):
        measured_concentration_colnames = [measured_concentration_colnames]
    if isinstance(concentration_in_feed, float) | isinstance(concentration_in_feed, int):
        concentration_in_feed = [concentration_in_feed]

    # validate input
    if len(measured_concentration_colnames) != len(concentration_in_feed):
        raise ValueError(
            "The number of species given in measured_concentration_colnames and "
            "concentration_in_feed does not match. Please check the input."
            "Remeber also add concentration of the species that are present in the "
            "feed, i.e. 0."
        )

    if isinstance(accumulated_feed_colname, Iterable) and not isinstance(accumulated_feed_colname, str):
        # if multiple feed streams are given, we validate that each list in the
        # the concentration_in_feed has the same length as the number of feeds
        for conc_lst in concentration_in_feed: 
            if not(isinstance(conc_lst, Iterable)):
                raise ValueError(
                    "If multiple feeds are given, the "
                    "concentration_in_feed argument should be a list of Iterables. "
                    f"{conc_lst} is not an Iterable. Please check the input."
                )

            if len(conc_lst) != len(accumulated_feed_colname):
                raise ValueError(
                    f"{len(accumulated_feed_colname)} feeds are given, but "
                    "one of the lists in concentration_in_feed has length "
                    f"{len(conc_lst)}. Please check the input."
                )


    out = pd.DataFrame()
    for species, conc in zip(
        measured_concentration_colnames, concentration_in_feed
    ):
        out[species + pseudo_col_postfix] = pseudobatch_transform(
            measured_concentration=df[species].to_numpy(),
            reactor_volume=df[reactor_volume_colname].to_numpy(),
            accumulated_feed=df[accumulated_feed_colname].to_numpy(),
            concentration_in_feed=conc,
            sample_volume=df[sample_volume_colname].to_numpy(),
        )
    # Copy the index from the original dataframe
    out.index = df.index
    return out


def convert_volumetric_rates_from_pseudo_to_real(
    pseudo_volumetric_rates: Union[pd.Series, NDArray],
    reactor_volume: Union[pd.Series, NDArray],
    sample_volume: Union[pd.Series, NDArray],
) -> Union[pd.Series, NDArray]:
    """
    Convert pseudo concentration to real concentration.

    Parameters
    ----------
    pseudo_concentration : Union[pd.Series, NDArray]
        Pseudo volumetric rates to be converted.
    reactor_volume : Union[pd.Series, NDArray]
        Reactor volume, this must be the BEFORE sampling volume.
    sample_volume : Union[pd.Series, NDArray]
        Sample volume.

    Returns
    -------
    Union[pd.Series, NDArray]
        Real volumetric rates.


    """
    adf = accumulated_dilution_factor(
        reactor_volume - sample_volume, sample_volume
    )
    return pseudo_volumetric_rates / adf


def hypothetical_concentration(
    metabolised_amount: np.ndarray,
    reactor_volume: np.ndarray,
    sample_volume: np.ndarray
) -> np.ndarray:
    """Preprocess the gaseous/volatile species data to prepare it for the pseudo batch
    transformation. The functions calculates the hypothetical liquid concentration
    of the gaseous species if the gas species did not evaporate. This is done by
    first calculating the hypothetical concentration of evaporated species including 
    losses due to sampling. Then the function adds the hypothetical concentration of
    evaporated species to the measured concentration of the gaseous species to get
    the preprocessed concentration data.

    Parameters
    ----------
    metabolised_amount : np.ndarray
        Net accumulated result of metabolism at each time point. Thus, if the species
        is consumed, the values are negative, if the species is produced, the values 
        are positive.
    reactor_volume : np.ndarray
        Reactor volume BEFORE sampling.
    sample_volume : np.ndarray
        Sample volume.
   
    Returns
    -------
    np.ndarray
        Hypothetical liquid concentration of gaseous species if the gas species
        did not evaporate.

    Notes
    -----
    For highly volatile species, such as CO2 and O2, the measured concentration
    in the liquid phase can be assumed to be zero. 
    """
    # Initialize an empty array to store the accumulated mass loss due to sampling.
    accumulated_amount_loss_due_to_sampling = np.zeros_like(metabolised_amount)

    # Calculate the accumulated mass loss due to sampling of the hypothetically 
    # dissolved gaseous molecules.
    for i in range(len(metabolised_amount)):
        if i == 0:
            # The first element of the array is calculated differently.
            accumulated_amount_loss_due_to_sampling[i] = metabolised_amount[i] * sample_volume[i] / reactor_volume[i]
        else:
            accumulated_amount_loss_due_to_sampling[i] = (
                (metabolised_amount[i] - accumulated_amount_loss_due_to_sampling[i-1]) * sample_volume[i] / reactor_volume[i]
                + accumulated_amount_loss_due_to_sampling[i-1]
            )

    # Calculate the sampling-adjusted gas concentration for each gaseous species.
    sampling_adjusted_concentration = (metabolised_amount - accumulated_amount_loss_due_to_sampling) / (reactor_volume - sample_volume)

    return sampling_adjusted_concentration

def metabolised_amount(
    off_gas_amount: NDArray,
    dissolved_amount_after_sampling: NDArray,
    inlet_gas_amount: NDArray,
    sampled_amount: NDArray,
    inlet_liquid_amount: Union[NDArray, None]= None,
):
    '''Calulated the amount of metabolised species at a given time.
    Essentially, solving the mass balance equation for the metabolised species.
    
    Parameters
    ----------
    off_gas_amount : NDArray
        Accumulated amount of the species exiting the reactor.
    dissolved_amount_after_sampling : NDArray
        Amount of the species in the reactor after sampling.
    inlet_gas_amount : NDArray
        Accumulated amount of the species entering the reactor through the gas inlet.
    sampled_amount : NDArray
        Accumulated amount of the species removed from the reactor through sampling.
    inlet_liquid_amount : NDArray, optional
        Accumulated amount of the species entering the reactor through the liquid inlet.
        If not given, it is assumed that no mass of the species entering the reactor through the
        liquid inlet.
    '''

    if inlet_liquid_amount is None:
        inlet_liquid_amount = np.zeros_like(off_gas_amount)

    return off_gas_amount + sampled_amount + dissolved_amount_after_sampling - inlet_gas_amount - inlet_liquid_amount 