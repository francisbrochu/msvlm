# __author__ = 'Alexandre Drouin'


class PreprocessorMixin:
    """
    A mixin class for the spectrum pre-processing algorithms.
    """

    def __init__(self):
        pass

    def fit(self, spectra_list):
        """
        Fit the pre-processing algorithm based on a training sample of spectra.

        Parameters
        ----------
        spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of training spectra.
        """
        pass

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
        raise NotImplementedError()

    def fit_transform(self, spectra_list):
        """
        Fit the pre-processing algorithm based on a training sample of spectra and transform the spectra based
        on the fitted parameters.

        Parameters
        ----------
        spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of spectra to transform.

        Returns
        -------
        transformed_spectra_list: array-like, type=Spectrum, shape=[n_spectra]
            The list of transformed spectra.
        """
        self.fit(spectra_list)
        return self.transform(spectra_list)
