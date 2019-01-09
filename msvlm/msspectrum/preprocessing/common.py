# __author__ = 'Alexandre Drouin'

import numpy as np
from ..base import PreprocessorMixin


class Pipeline(PreprocessorMixin):
    """
    A pipeline class for applying multiple pre-processing algorithms sequentially.
    """
    def __init__(self, preprocessors):
        """
        Constructor

        Parameters
        ----------
        preprocessors: array-like, type=any_pre-processing_class, shape=[n_pre-processing_steps]
            A list of pre-processing objects to apply sequentially.
        """
        self.preprocessors = preprocessors

    def fit(self, spectra_list):
        """
        Fit the pipeline based on a training sample of spectra.

        Each preprocessor is fitted and used to transform the training sample, which is subsequently used to fit the
        next preprocessor.

        Parameters
        ----------
        spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of training spectra.
        """
        spectra_list = np.array(spectra_list)
        for preprocessor in self.preprocessors:
            spectra_list = preprocessor.fit_transform(spectra_list)

    def transform(self, spectra_list):
        """
        Transform a list of spectra based on the fitted parameters.

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
        for preprocessor in self.preprocessors:
            spectra_list = preprocessor.transform(spectra_list)

        return spectra_list
