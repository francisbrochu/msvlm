# -*- coding: utf-8 -*-

import numpy as np
from bisect import bisect_left
from .msspectrum.preprocessing.discrete import ThresholdedPeakFiltering, MassRangeSelection
from .msspectrum.spectrum import Spectrum
from .msspectrum.utils import binary_search_for_right_range, binary_search_for_left_range
from copy import deepcopy
import msvlm.msAlign.msAlign as ms


def take_closest(my_list, my_number, lo=0):
    # http://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    """
    Assumes my_list is sorted. Returns closest value to my_number.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(my_list, my_number, lo=lo)
    if pos == 0:
        return my_list[0], pos
    if pos == len(my_list):
        return my_list[-1], pos - 1
    before = my_list[pos - 1]
    after = my_list[pos]
    if after - my_number < my_number - before:
        return after, pos
    else:
        return before, pos-1


def is_window_vlm(peaks, window_start_idx, window_end_idx, w):
    """
    TODO

    """
    start_mass = peaks[window_start_idx]
    end_mass = peaks[window_end_idx]
    center_mass = np.mean(peaks[window_start_idx:window_end_idx + 1])
    window_lower_bound = center_mass * (1 - w)
    window_upper_bound = center_mass * (1 + w)
    if(window_lower_bound > start_mass) or (window_upper_bound < end_mass):
        return False
    # start mass and end mass are within window

    # check if peak before and after the indexes are within the window, if yes then not a VLM
    mass_lower = 0.
    if window_start_idx != 0:
        mass_lower = peaks[window_start_idx - 1]

    mass_over = 1000000.
    if window_end_idx != (len(peaks) - 1):
            mass_over = peaks[window_end_idx + 1]

    if (mass_lower >= window_lower_bound) or (mass_over <= window_upper_bound):
        return False

    return True


class VirtualLockMassCorrector(object):
    def __init__(self, window_size, minimum_peak_intensity, max_skipped_points=None):
        """
        Initiate a VirtualLockMassCorrector object.
        :param window_size: The distance from left to right in ppm
        :param minimum_peak_intensity: Minimum peak intensity to be considered by the algorithm
        :param max_skipped_points: Maximum number of points that can be skipped during the transform step. None=any.
        """
        self._is_fitted = False
        self.window_size = window_size
        self.window_size_ppm = 1.0 * window_size / 10**6
        self.minimum_peak_intensity = minimum_peak_intensity
        self._vlm_mz = None
        self.max_skipped_points = max_skipped_points
        self._used_vlm = None
        self._mz_precision = 4

    def _preprocess_spectra(self, spectra):
        return ThresholdedPeakFiltering(threshold=self.minimum_peak_intensity,
                                        remove_mz_values=True).fit_transform(spectra)

    def _find_vlm_peaks(self, spectra):
        spectra = self._preprocess_spectra(spectra)
        peaks_list = []
        for s in spectra:
            mzs = [mz for mz in s.mz_values]
            peaks_list.append(mzs)

        vlm_mz_values = ms.find_vlm(peaks_list, self.window_size * 10**-6)
        return vlm_mz_values

    def _apply_correction(self, spectrum):
        """
        Apply the VLM to a spectrum
        :param spectrum: A pymspec spectrum to correct
        :return: A corrected pymspec spectrum.
        """
        if len(self._vlm_mz) <= 3:
            raise ValueError("There must be at least 3 points to use virtual lock-mass")

        spectrum_copy = deepcopy(spectrum)

        # Find the corresponding points
        found_vlm, observed_mz = self._find_vlock_mass_in_spectra(spectrum_copy)

        # Cut peaks of m/z lower than first VLM and greater than last VLM
        selector = MassRangeSelection(lower_range=observed_mz[0], upper_range=observed_mz[-1])
        spectrum_copy = selector.fit_transform([spectrum])[0]
        spect_mz = spectrum_copy.mz_values

        corrected_mzs = []
        idx = 0

        while idx < (len(found_vlm) - 1):
            corrected_mzs += self._correct_points_between(spect_mz, observed_mz[idx], observed_mz[idx + 1],
                                                          found_vlm[idx], found_vlm[idx + 1])
            idx += 1

        corrected_mzs += [found_vlm[-1]]  # need to add final VLM

        new_spectrum = Spectrum(mz_values=np.asarray(corrected_mzs),
                                intensity_values=spectrum_copy.intensity_values,
                                mz_precision=spectrum_copy.mz_precision,
                                metadata=spectrum_copy.metadata)
        return new_spectrum

    def _find_vlock_mass_in_spectra(self, spectrum):
        """
        Search each vlm in a spectrum and return the list of vlm found and their correspondance in the spectrum.
        :param spectrum: A pymspec spectrum
        :return: two lists: The vlm found in the spectrum and their correspondance
        """
        preprocessed_spect = self._preprocess_spectra([spectrum])[0]

        observed_mz = []
        vlm_found = []
        number_skipped_points = 0
        for vlm in self._vlm_mz:
            last_index = 0
            try:
                best_match, position = take_closest(preprocessed_spect.mz_values, vlm,
                                                    lo=last_index)
                # Check if the vlm is in the window
                mz_difference = abs(best_match - vlm)
                # self.window_size is from center to side, not side to side.
                if mz_difference > vlm * self.window_size_ppm:
                    raise ValueError("A VLM was not found in the appropriate window")
                last_index = position
                observed_mz.append(best_match)
                vlm_found.append(vlm)
            except ValueError as error:
                if self.max_skipped_points is None:
                    pass  # If none, any VLM can not be found
                else:
                    number_skipped_points += 1
                    if number_skipped_points > self.max_skipped_points:
                        raise error
        observed_mz = np.array(observed_mz)
        vlm_found = np.array(vlm_found)
        np.around(observed_mz, decimals=self._mz_precision)
        return vlm_found, observed_mz

    def _correct_points_between(self, mzs, observed_mz1, observed_mz2, ref_mz1, ref_mz2):
        """
        :param observed_mz1: an observed mz of a virtual lock mass (smaller than the 2nd)
        :param observed_mz2: an observed mz of a virtual lock mass (greater than the first)
        :param ref_mz1: a reference mz of a virtual lock mass (smaller than the 2nd)
        :param ref_mz2: a reference mz of a virtual lock mass (greater than the first)
        :return: the corrected mz values from the spectrum that are between mz1 and mz2
        """
        if observed_mz1 <= 0 or observed_mz2 <= 0 or ref_mz1 <= 0 or ref_mz2 <= 0:
            raise ValueError("Mz and ratios cannot be null or negative")
        function = self._create_correction_function(observed_mz1, observed_mz2, ref_mz1, ref_mz2)

        right = binary_search_for_right_range(mzs, observed_mz2)
        left = binary_search_for_left_range(mzs, observed_mz1)
        mz_to_be_corrected = mzs[left:right]
        corrected_mz = [function(mz) for mz in mz_to_be_corrected]
        return corrected_mz

    def _create_correction_function(self, mz1, mz2, ref1, ref2):
        """
        Create the y = m*x + b function for the 2 points in parameter
        :param mz1: lowest mz
        :param mz2: highest mz
        :param ref1: reference at mz1
        :param ref2: reference ratio at mz2
        :return: a numpy function that can correct the values between mz1 and mz2
        :raises: ValueError is an mz or ratio is <= 0
        """
        if mz1 <= 0 or mz2 <= 0 or ref1 <= 0 or ref2 <= 0:
            raise ValueError("Mz cannot be null or negative")
        if mz1 > mz2:
            raise ValueError("mz2 must be greater than mz1")
        m = (ref2 - ref1) / (mz2 - mz1)
        b = ref2 - (m * mz2)
        function = lambda x: m * x + b

        return function

    def optimize_window_size(self, spectrum,
                             possible_values=(5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)):
        estimation_results = []
        for p in possible_values:
            vlm_estimator = VirtualLockMassCorrector(p, self.minimum_peak_intensity)
            vlm_estimator.fit(spectrum)
            estimation_results.append(len(vlm_estimator._vlm_mz))
            if len(estimation_results) > 3:
                if estimation_results[-1] < estimation_results[-3] > estimation_results[-2]:
                        break
        best_p = max(estimation_results)
        index = estimation_results.index(best_p)
        self.window_size = possible_values[index]
        self.window_size_ppm = 1.0 * self.window_size / 10**6

    def fit(self, spectra):
        """
        TODO

        """
        self._mz_precision = spectra[0].mz_precision
        self._vlm_mz = self._find_vlm_peaks(spectra)
        self._is_fitted = True

    def transform(self, spectra):
        """
        TODO

        """
        if self._vlm_mz is None:
            raise RuntimeError("The VLM corrector must be fitted before applying a correction.")
        return np.asarray([self._apply_correction(spectrum) for spectrum in spectra])

    def _choose_vlm_points(self, n_skip=1):
        used_vlm = np.zeros(shape=len(self._vlm_mz), dtype=bool)

        for i, vlm in enumerate(self._vlm_mz):
            if (i + 1) % (1 + n_skip) == 0:
                used_vlm[i] = False
            else:
                used_vlm[i] = True

        used_vlm[0] = True
        used_vlm[-1] = True

        return used_vlm

    def _tag_peaks(self, spectrum, observed_vlms):
        vlm_tag = np.zeros(shape=len(spectrum.mz_values), dtype=bool)

        vlm_idx = 0

        i = 0
        while (i < len(spectrum.mz_values)) and (vlm_idx < len(observed_vlms)):

            if spectrum.mz_values[i] == observed_vlms[vlm_idx]:
                vlm_tag[i] = True
                vlm_idx += 1

            i += 1

        return vlm_tag

    def _apply_loo_correction(self, spectrum):

        num_vlm_used = np.sum(self._used_vlm)
        if num_vlm_used < 3:
            raise ValueError("There must be at least 3 points to use virtual lock-mass")

        spectrum_copy = deepcopy(spectrum)

        # Find the corresponding points
        found_vlm, observed_mz = self._find_vlock_mass_in_spectra(spectrum_copy)
        # since we fit and _apply_optimize_alignment on the same set of spectra, all VLMs will be found
        # mark them in the spectrum
        spectrum_vlm_tags = self._tag_peaks(spectrum_copy, observed_mz)
        # Use only selected VLMs
        use_found_vlm = found_vlm[self._used_vlm]
        use_observed_mz = observed_mz[self._used_vlm]

        vlm_spectrum = Spectrum(mz_values=spectrum_copy.mz_values[spectrum_vlm_tags],
                                intensity_values=spectrum_copy.intensity_values[spectrum_vlm_tags],
                                mz_precision=spectrum_copy.mz_precision,
                                metadata=spectrum_copy.metadata)

        # Cut peaks of m/z lower than first VLM and greater than last VLM
        selector = MassRangeSelection(lower_range=use_observed_mz[0], upper_range=use_observed_mz[-1])
        vlm_spectrum = selector.fit_transform([vlm_spectrum])[0]
        vlm_spect_mz = vlm_spectrum.mz_values

        corrected_vlm_mzs = []
        idx = 0

        while idx < (len(use_found_vlm) - 1):
            corrected_vlm_mzs += self._correct_points_between(vlm_spect_mz, use_observed_mz[idx],
                                                              use_observed_mz[idx + 1], use_found_vlm[idx],
                                                              use_found_vlm[idx + 1])
            idx += 1

        corrected_vlm_mzs += [use_found_vlm[-1]]  # need to add final VLM

        new_vlm_spectrum = Spectrum(mz_values=np.asarray(corrected_vlm_mzs),
                                    intensity_values=vlm_spectrum.intensity_values,
                                    mz_precision=vlm_spectrum.mz_precision,
                                    metadata=vlm_spectrum.metadata)

        return new_vlm_spectrum

    def leave_one_out_optimize(self, spectra, proportion=0.95):
        if not self._is_fitted:
            raise RuntimeError("The VLM corrector must be fitted before applying a correction.")

        vlm_max_dist = []
        for i in range(len(self._vlm_mz) - 2):
            idx = i + 1

            vlm_mzs = []
            for spectrum in spectra:
                found_vlm, observed_mz = self._find_vlock_mass_in_spectra(spectrum)
                vlm_mzs.append(observed_mz)
            vlm_mzs = np.asarray(vlm_mzs)

            # each vlm "spectrum" should be of same length, the number of VLMs, as this is applied on a training set

            mzs = []
            for j in range(len(vlm_mzs)):
                corr_function = self._create_correction_function(vlm_mzs[j][idx - 1], vlm_mzs[j][idx + 1],
                                                                 self._vlm_mz[idx - 1], self._vlm_mz[idx + 1])
                mzs.append(corr_function(vlm_mzs[j][idx]))
            mzs = np.asarray(mzs)
            comp_mz = np.mean(mzs)

            dist = np.abs(mzs - comp_mz)
            ppm_dist = dist * 10**6 / comp_mz
            vlm_max_dist.append(np.max(ppm_dist))

        sorted_distances = np.sort(vlm_max_dist)
        idx = int(proportion * len(sorted_distances)) - 1
        chosen_dist = sorted_distances[idx]

        return chosen_dist, vlm_max_dist
