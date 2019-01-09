# -*- coding: utf-8 -*-

import numpy as np
from .msspectrum.spectrum import Spectrum
from .msspectrum.utils import take_closest, binary_search_mz_values
from heapq import heappop, heappush, heapreplace
from msvlm.msAlign.msAlign import find_alpt


class MassSpectraAligner:

    def __init__(self, window_size=10):
        self.window_size = window_size
        self.reference_mz = []

    def fit(self, spectra):
        self._train(spectra)

    def _train(self, spectra):
        """
        Fill the reference_mz attribute with possible m/z values.
        :param spectra: A set of spectrum object.
        :return: Nothing
        """
        peaks_list = []
        for s in spectra:
            mzs = [mz for mz in s.mz_values]
            peaks_list.append(mzs)

        self.reference_mz = find_alpt(peaks_list, self.window_size * 10**-6)

    def transform(self, spectra):
        new_spectra = []
        for s in spectra:
            new_spectra.append(self._apply(s))
        return np.asarray(new_spectra)

    def transform_test(self, spectra):
        new_spectra = []
        for s in spectra:
            new_spectra.append(self._apply_test(s))
        return np.asarray(new_spectra)

    def _apply(self, spec):
        # Find closest point that is not outside possible window
        # If point: change mz
        # Else: keep or discard m/z?
        aligned_mz = []
        aligned_int = []
        nf_mz = []
        nf_int = []
        for i, mz in enumerate(spec.mz_values):
            try:
                possible_matches = binary_search_mz_values(self.reference_mz, mz,
                                                           float(self.window_size))
            except ValueError:
                nf_mz.append(mz)
                nf_int.append(spec.intensity_values[i])
                continue

            if len(possible_matches) > 1:
                possible_matches = [take_closest(possible_matches, mz)]

            if len(possible_matches) == 1:
                aligned_mz.append(possible_matches[0])
                aligned_int.append(spec.intensity_values[i])
            else:
                aligned_mz.append(mz)
                aligned_int.append(spec.intensity_values[i])
                nf_mz.append(mz)
                nf_int.append(spec.intensity_values[i])

        return Spectrum(np.asarray(aligned_mz), np.asarray(aligned_int),
                        spec.mz_precision, spec.metadata)

    def _apply_test(self, spec):
        # Find closest point that is not outside possible window
        # If point: change mz
        # Else: keep or discard m/z?
        aligned_mz = []
        aligned_int = []
        for i, mz in enumerate(spec.mz_values):
            try:
                possible_matches = binary_search_mz_values(self.reference_mz, mz,
                                                           float(self.window_size))
            except ValueError:
                continue

            if len(possible_matches) > 1:
                possible_matches = [take_closest(possible_matches, mz)]

            if len(possible_matches) == 1:
                aligned_mz.append(possible_matches[0])
                aligned_int.append(spec.intensity_values[i])

        return Spectrum(np.asarray(aligned_mz), np.asarray(aligned_int),
                        spec.mz_precision, spec.metadata)

    def optimize_window(self, spectra, vlm_spectra,
                        window=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 18, 20, 22, 25, 28, 30)):

        res = np.zeros(shape=(len(window)))
        n_align = np.zeros(shape=(len(window)))
        align_populations = []
        refs = []

        # fit with each window size
        for i, w in enumerate(window):
            opt_aligner = MassSpectraAligner(window_size=w)
            opt_aligner.fit(spectra)
            res[i], tmp1, tmp2 = opt_aligner._evaluate_alignment(vlm_spectra)
            n_align[i] = len(opt_aligner.reference_mz)
            align_populations.append((tmp1, tmp2))
            refs.append(opt_aligner.reference_mz)

        align_populations = np.asarray(align_populations)

        max_score = 0
        best_window_size = 1
        for i, p in enumerate(align_populations):
            score = 0
            for el in p[1]:
                if el == len(vlm_spectra):
                    score += 1

            if score > max_score:
                max_score = score
                best_window_size = window[i]

        return best_window_size, align_populations, refs

    def _evaluate_alignment(self, spectra):
        a_idx = 0
        a_mz = self.reference_mz[a_idx]
        a_lb = a_mz - (a_mz * self.window_size / 1000000.)
        a_ub = a_mz + (a_mz * self.window_size / 1000000.)
        s_rep = {}

        h = []
        spectra_idx = np.zeros(shape=len(spectra), dtype=int)
        alignment_status = np.zeros(len(self.reference_mz), dtype=bool)
        a_pop = np.zeros(len(self.reference_mz), dtype=int)
        spectra_max = np.zeros(shape=len(spectra), dtype=int)

        for i in range(len(spectra)):
            heappush(h, (spectra[i].mz_values[spectra_idx[i]], i))
            spectra_idx[i] += 1
            spectra_max[i] = len(spectra[i].mz_values)

        while (len(h) != 0) and (a_idx < len(self.reference_mz)):
            next_s = h[0][1]

            if spectra_idx[next_s] < spectra_max[next_s]:
                mz, s = heapreplace(h, (spectra[next_s].mz_values[spectra_idx[next_s]], next_s))  # advance in s
                spectra_idx[next_s] += 1
            else:  # last peak of a spectra next_s, remove it from heap
                mz, s = heappop(h)

            if mz < a_lb:  # peak before current alignment point,
                pass
            elif mz > a_ub:  # peak after current alignment point
                s_rep = {}
                # advance until the current alignment point is not passed by current peak
                while (mz > a_ub) and (a_idx < len(self.reference_mz) - 1):
                    a_idx += 1
                    a_mz = self.reference_mz[a_idx]
                    a_lb = a_mz - (a_mz * self.window_size / 1000000.)
                    a_ub = a_mz + (a_mz * self.window_size / 1000000.)

            # if peak is within current alignment point
            if (mz >= a_lb) and (mz <= a_ub):
                if len(s_rep) == 0:  # first peak in an alignment point window
                    alignment_status[a_idx] = True
                    s_rep[s] = True
                    a_pop[a_idx] += 1
                else:
                    # check if spectrum already in alignment window
                    if s in s_rep:
                        alignment_status[a_idx] = False
                        a_pop[a_idx] += 1
                    else:
                        s_rep[s] = True
                        a_pop[a_idx] += 1

        return np.sum(alignment_status), alignment_status, a_pop
