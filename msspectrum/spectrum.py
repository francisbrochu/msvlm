# __author__ = 'Alexandre Drouin'
from .utils import _is_mz_precision_equal
import numpy as np
from copy import deepcopy


class Spectrum(object):
    def __init__(self, mz_values, intensity_values, mz_precision=3, metadata=None):
        self._peaks_mz = np.array([])
        self._peaks_intensity = np.array([])
        self._peaks = {}
        self.metadata = metadata
        self._mz_precision = mz_precision  # in decimals e.g.: mz_precision=3 => 5.342

        if len(mz_values) != len(intensity_values):
            raise ValueError("The number of mz values must be equal to the number of intensity values.")

        self.set_peaks(mz_values, intensity_values)

    def peaks(self):
        """
        Note: Peaks are not necessarily sorted here because of dict
        """
        return self._peaks

    @property
    def mz_values(self):
        """
        Note: Returned values are always sorted
        """
        return self._peaks_mz

    @property
    def mz_precision(self):
        return self._mz_precision

    @mz_precision.setter
    def mz_precision(self, new_precision):
        self._mz_precision = new_precision
        self.set_peaks(self.mz_values, self.intensity_values)

    @property
    def intensity_values(self):
        return self._peaks_intensity

    def intensity_at(self, mz):
        mz = round(mz, self._mz_precision)
        try:
            intensity = self._peaks[mz]
        except AttributeError:
            intensity = 0.0
        return intensity

    def set_peaks(self, mz_values, intensity_values):
        # XXX: This function must create a copy of mz_values and intensity_values to prevent the modification of
        # referenced arrays. This is assumed by other functions. Be careful!

        # Sort the peaks by mz
        sort_mz = np.argsort(mz_values)
        mz_values = np.asarray(mz_values)[sort_mz]
        intensity_values = np.asarray(intensity_values)[sort_mz]

        # Round the mz values based on the mz precision
        mz_values = np.asarray(np.round(mz_values, self._mz_precision), dtype=np.float)

        # Contiguous mz values might now be equivalent. Combine their intensity values by taking the sum.
        unique_mz = np.unique(mz_values)
        unique_mz_intensities = np.zeros(unique_mz.shape)

        # Note: This assumes that mz_values and unique_mz are sorted
        if len(mz_values) != len(unique_mz):
            acc = 0
            current_unique_mz_idx = 0
            current_unique_mz = unique_mz[0]
            for i, mz in enumerate(mz_values):
                if mz != current_unique_mz:
                    unique_mz_intensities[current_unique_mz_idx] = acc  # Flush the accumulator
                    acc = 0  # Reset the accumulator
                    current_unique_mz_idx += 1  # Go to the next unique mz value
                    current_unique_mz = unique_mz[current_unique_mz_idx]  # Get the unique mz value
                acc += intensity_values[i]  # Increment the accumulator
            unique_mz_intensities[current_unique_mz_idx] = acc  # Flush the accumulator
        else:
            unique_mz_intensities = intensity_values

        self._peaks_mz = unique_mz
        self._peaks_mz.flags.writeable = False

        self._peaks_intensity.flags.writeable = False
        self._peaks_intensity = unique_mz_intensities

        self._peaks = dict([(round(self._peaks_mz[i], self._mz_precision), self._peaks_intensity[i]) for i in
                            range(len(self._peaks_mz))])

        self._check_peaks_integrity()

    def copy(self):
        return copy_spectrum(self)

    def __iter__(self):
        """
        Returns an iterator on the peaks of the spectrum.

        Returns:
        --------
        peak_iterator: iterator
            An iterator that yields tuples of (mz, int) for each peaks in the spectrum.
        """
        return zip(self._peaks_mz, self._peaks_intensity)

    def __len__(self):
        return self._peaks_mz.shape[0]

    def _check_peaks_integrity(self):
        if not len(self._peaks_mz) == len(self._peaks_intensity):
            raise ValueError("The number of mz values must be equal to the number of intensity values.")
        if not all(self._peaks_mz[i] <= self._peaks_mz[i + 1] for i in range(len(self._peaks_mz) - 1)):
            raise ValueError("Mz values must be sorted.")
        if len(np.unique(self._peaks_mz)) != len(self._peaks_mz):
            raise ValueError("Mz value list contains duplicate values.")


def copy_spectrum(spectrum):
    """
    Copies a spectrum.

    Parameters:
    -----------
    spectrum: Spectrum
        The spectrum to copy

    Note:
    -----
    * This ensures that the metadata is deepcopied
    """
    metadata = deepcopy(spectrum.metadata)
    # XXX: The mz_values and intensity_values are copied in the constructor. No need to copy here.
    return Spectrum(mz_values=spectrum.mz_values, intensity_values=spectrum.intensity_values,
                    mz_precision=int(spectrum.mz_precision), metadata=metadata)


def copy_spectrum_with_new_intensities(spectrum, new_intensity_values):
    """
    Copies a spectrum and replaces its intensity values.

    Parameters:
    -----------
    spectrum: Spectrum
        The spectrum to copy
    new_intensity_values: array_like, dtype=float, shape=n_peaks
        The new intensity values

    Note:
    -----
    * This is more efficient than deepcopying the spectrum and modifying its intensity values.
    * This ensures that the metadata is deepcopied
    """
    metadata = deepcopy(spectrum.metadata)
    # XXX: The mz_values and intensity_values are copied in the constructor. No need to copy here.
    return Spectrum(mz_values=spectrum.mz_values, intensity_values=new_intensity_values,
                    mz_precision=int(spectrum.mz_precision), metadata=metadata)


def copy_spectrum_with_new_mz_and_intensities(spectrum, new_mz_values, new_intensity_values):
    """
    Copies a spectrum and replaces its mz and intensity values.

    Parameters:
    -----------
    spectrum: Spectrum
        The spectrum to copy
    new_mz_values: array_like, dtype=float, shape=n_peaks
        The new mz values
    new_intensity_values: array_like, dtype=float, shape=n_peaks
        The new intensity values

    Note:
    -----
    * This is more efficient than deepcopying the spectrum and modifying its mz and intensity values.
    * This ensures that the metadata is deepcopied
    """
    metadata = deepcopy(spectrum.metadata)
    # XXX: The mz_values and intensity_values are copied in the constructor. No need to copy here.
    return Spectrum(mz_values=new_mz_values, intensity_values=new_intensity_values,
                    mz_precision=int(spectrum.mz_precision), metadata=metadata)


def unify_mz(spectra):
    """
    Unifies the m/z values for a list of spectra

    Parameters:
    -----------
    spectra: list of Spectrum
        A list of spectra.

    Note:
    -----
    * The operation is performed in-place
    """
    if not _is_mz_precision_equal(spectra[0].mz_precision, spectra):
        raise ValueError("The m/z precision of the spectra must be equal in order to unify the m/z values.")

    mz_values = np.array(sorted(set(mz for spectra in spectra for mz in spectra.mz_values)))

    for spectrum in spectra:
        spectrum.set_peaks(mz_values=mz_values,
                           intensity_values=np.array([spectrum.intensity_at(mz) for mz in mz_values]))


def unify_precision(spectra, new_precision):
    """
    Unifies the m/z precision for a list of spectra

    Parameters:
    -----------
    spectra: list of Spectrum
        A list of spectra.

    new_precision: int
        The number of decimals for the new precision

    Note:
    -----
    * The operation is performed in-place
    * If multiple m/z values are equal after adjusting the precision, their intensity values are summed.
    """
    for spectrum in spectra:
        spectrum.mz_precision = new_precision
