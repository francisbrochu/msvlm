# __author__ = 'Alexandre Drouin'

from ..base import PreprocessorMixin
from ..spectrum import copy_spectrum_with_new_intensities, copy_spectrum_with_new_mz_and_intensities
from ..utils import binary_search_for_left_range, binary_search_for_right_range
import copy
import numpy as np


class ThresholdedPeakFiltering(PreprocessorMixin):
    """
    A pre-processor for removing the peaks that are less intense than a given threshold.
    """
    def __init__(self, threshold=1.0, remove_mz_values=True):
        """
        Constructor.

        Parameters
        ----------
        threshold: float, default=1.0
                   The intensity threshold. All peaks that have an intensity value less or equal to this threshold will
                   be discarded.

        remove_mz_values : bool, default=True
                   Specifies if the m/z values where the intensity was below the threshold should be removed.
        """
        self.threshold = threshold
        self.remove_mz_values = remove_mz_values

    def transform(self, spectra_list):
        """
        Filter peaks for a list of spectra.

        Parameters
        ----------
        spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of spectra to transform.

        Returns
        -------
        transformed_spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of transformed spectra.
        """

        spectra_list = np.array(spectra_list)
        for i, spectrum in enumerate(spectra_list):
            if not self.remove_mz_values:
                intensity_values = copy.deepcopy(spectra_list[i].intensity_values)
                intensity_values[intensity_values <= self.threshold] = 0.0
                spectra_list[i] = copy_spectrum_with_new_intensities(spectrum, intensity_values)
            else:
                keep_mask = spectra_list[i].intensity_values > self.threshold
                spectra_list[i] = copy_spectrum_with_new_mz_and_intensities(spectrum,
                                                                            spectra_list[i].mz_values[keep_mask],
                                                                            spectra_list[i].intensity_values[keep_mask])
        return spectra_list


class MassRangeSelection(PreprocessorMixin):
    """
    A pre-processor for removing the peaks that are outside a given window on the m/z range
    """
    def __init__(self, lower_range=50.0, upper_range =2000.0):
        """
        Constructor.

        Parameters
        ----------
        lower_range: float, default=50.0
                   The intensity threshold. All peaks that have an intensity value less or equal to this threshold will
                   be discarded.

        upper_range : float, default=2000.0
                   Specifies if the m/z values where the intensity was below the threshold should be removed.
        """
        self.lower_range = lower_range
        self.upper_range = upper_range

    def transform(self, spectra_list):
        """
        Filter peaks for a list of spectra.

        Parameters
        ----------
        spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of spectra to transform.

        Returns
        -------
        transformed_spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of transformed spectra.
        """
        spectra_list = np.array(spectra_list)
        for i, spectrum in enumerate(spectra_list):
            mz_values = spectra_list[i].mz_values
            intensity_values = spectra_list[i].intensity_values

            # mz_values is sorted and intensity values in the same order.
            # find the indices of the peaks are within the range and slice the arrays.
            lower_idx = binary_search_for_left_range(mz_values, self.lower_range)
            upper_idx = binary_search_for_right_range(mz_values, self.upper_range) + 1

            spectra_list[i] = copy_spectrum_with_new_mz_and_intensities(spectrum, mz_values[lower_idx:upper_idx],
                                                                        intensity_values[lower_idx:upper_idx])

        return spectra_list
