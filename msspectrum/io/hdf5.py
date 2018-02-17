from __future__ import print_function, division, absolute_import, unicode_literals
__author__ = 'Alexandre Drouin'


import h5py as h
import json

from numpy import float

from ..spectrum import Spectrum
from ..utils import _is_mz_equal, _is_mz_precision_equal, _is_metadata_type_dict

def hdf5_save(file_name, spectra, compression_type="gzip", compression_opts=4):
    """
    Saves a list of spectra to a HDF5 file.

    Parameters:
    -----------
    file_name: str
        The path to the output file.

    spectra: list of Spectrum
        The list of spectra to save.

    compression_type: str, default="gzip"
        The type of compression to use for compressing the datasets. Refer to H5py's documentation.

    compression_opts: variable, default=4
        The compression options for the specified compression type. Refer to H5py's documentation. In the case of gzip
        compression, this specifies the level of compression used from 0 to 9.

    Note:
    -----
    * The spectra must have identical m/z values.
    * The spectra must have the same m/z precision.
    """
    if not _is_mz_equal(spectra[0].mz_values, spectra):
        raise ValueError("The spectra m/z values must be identicial to save to HDF5.")

    if not _is_mz_precision_equal(spectra[0].mz_precision, spectra):
        raise ValueError("The spectra m/z precisions must be equal to save to HDF5.")

    if not _is_metadata_type_dict(spectra):
        raise ValueError("The spectra metadata must be of type dict or None to save to HDF5.")

    file = h.File(file_name, "w")

    n_spectra = len(spectra)
    n_mz_values = len(spectra[0])

    file.create_dataset("precision", data=spectra[0].mz_precision)

    file.create_dataset("mz", data=spectra[0].mz_values, compression=compression_type,
                        compression_opts=compression_opts)

    spectra_intensity_dataset = file.create_dataset("intensity",
                                                    shape=(n_spectra, n_mz_values),
                                                    dtype=float,
                                                    chunks=(1, n_mz_values),
                                                    compression=compression_type,
                                                    compression_opts=compression_opts)

    for i, spectrum in enumerate(spectra):
        spectra_intensity_dataset[i] = spectrum.intensity_values

    dt = h.special_dtype(vlen=str)
    spectra_metadata_dataset = file.create_dataset("metadata", shape=(n_spectra,), dtype=dt)

    for i, spectrum in enumerate(spectra):
        spectra_metadata_dataset[i] = json.dumps(spectrum.metadata)

    file.close()


def hdf5_load(file_name, metadata=True):
    """
    Loads spectra from a HDF5 file.

    Parameters:
    -----------
    file_name: str
        The path to the file to load.

    metadata: boolean
        Defaults to True. Boolean to check if we load the metadata along with the spectrum data.

    Returns:
    -------
    spectra: list of Spectrum
        The list of spectra extracted form the file.
    """
    file = h.File(file_name, "r")

    mz_precision = file["precision"][...]
    mz_values = file["mz"][...]
    spectra_intensity_dataset = file["intensity"]

    if metadata and "metadata" in file:
        spectra_metadata_dataset = file['metadata']
    else:
        spectra_metadata_dataset = [None] * spectra_intensity_dataset.shape[0]

    spectra = []
    for spectrum_intensity_values, spectrum_metadata in zip(spectra_intensity_dataset, spectra_metadata_dataset):
        try:
            spectrum_metadata = json.loads(spectrum_metadata.decode("utf-8")) if not spectrum_metadata is None else None
        except AttributeError:
            spectrum_metadata = json.loads(spectrum_metadata) if not spectrum_metadata is None else None
        spectra.append(Spectrum(mz_values=mz_values, intensity_values=spectrum_intensity_values,
                                mz_precision=mz_precision, metadata=spectrum_metadata))
    file.close()

    return spectra
