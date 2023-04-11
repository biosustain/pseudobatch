import pandas as pd
import numpy as np
from typing import Union, Iterable
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
        Iterable[float], Iterable[NDArray[np.float64]]
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
        out[species + pseudo_col_postfix] = pseudobatch_transform(
            measured_concentration=df[species].to_numpy(),
            reactor_volume=df[reactor_volume_colname].to_numpy(),
            accumulated_feed=df[accumulated_feed_colname].to_numpy(),
            concentration_in_feed=conc,
            sample_volume=df[sample_volume_colname].to_numpy(),
        )
    return out

def reverse_pseudobatch_transform(
        pseudo_concentration: NDArray, 
        reactor_volume : NDArray,
        accumulated_feed: NDArray, 
        concentration_in_feed: Union[NDArray, float], 
        sample_volume:NDArray
) -> NDArray:
    '''Does the reverse of the pseudobatch transform, i.e. transforms pseudo 
    concentrations to real concentrations. The reverse transform is simply a
    rearrangement of the pseudobatch transformation formula were the concentration
    is isolated rather than the pseudo concentration. 

    Parameters
    ----------
    pseudo_concentration: NDArray
        a np.array with the pseudo concentrations of the species that should be transformed, e.g. biomass or compound.
        reactor_volume : NDArray

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
        A NDArray with the reverse transformed concetrations.

    Notes
    -----
    The reverse pseudo batch transformation formula is:

    .. math::

        \frac{G^{\star}_k + \sum_{i=1}^{k-1}ADF_i\frac{C_{feed}^{glucose}\cdot 
        F_{i}}{V_i}}{ADF_k} = C^{glucose}_k

    '''
    after_sample_reactor_volume = reactor_volume - sample_volume

    for i in [
        pseudo_concentration,
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
        return pseudo_concentration + np.cumsum(
            fed_species_term(
                accumulated_feed, concentration_in_feed, reactor_volume
            ) / adf
        )

    # Iterate over the columns of the accumulated feed to calculate
    # fed_species_term for each feed. Then calculate the row-wise sum of the
    # fed_species_term.

    fed_species_term_collection = np.zeros(np.shape(accumulated_feed)[0])
    for i in range(np.shape(accumulated_feed)[1]):
        fed_species_term_collection += fed_species_term(
            accumulated_feed[:, i], concentration_in_feed[i], reactor_volume
        )

    return pseudo_concentration + np.cumsum(fed_species_term_collection) / adf