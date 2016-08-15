# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio
import numpy as np


def _most_conserved(frequencies, sequence_dtype, gap_mode='ignore'):
    if gap_mode != 'ignore':
        raise ValueError('Unknown gap_mode: %s. ignore is currently the only '
                         'supported gap_mode.' % gap_mode)
    result = []
    for frequency_vector in frequencies:
        for gap in sequence_dtype.gap_chars:
            if gap in frequency_vector:
                del frequency_vector[gap]
        frequency_values = list(frequency_vector.values())
        if len(frequency_values) == 0:
            result.append(0.0)
        else:
            result.append(np.max(frequency_values) / np.sum(frequency_values))
    return result


def _compute_conservation_mask(frequencies, sequence_dtype, min_conservation):
    mask = [c >= min_conservation
            for c in _most_conserved(frequencies, sequence_dtype)]
    return np.array(mask)


def _compute_gap_mask(frequencies, sequence_dtype, max_gap_frequency):
    gap_frequencies = []
    num_sequences = np.sum(list(frequencies[0].values()))
    for f in frequencies:
        gap_frequency = np.sum([f.get(gc, 0.0)
                                for gc in sequence_dtype.gap_chars])
        gap_frequencies.append(gap_frequency/num_sequences)
    mask = [f <= max_gap_frequency for f in gap_frequencies]
    return np.array(mask)


def _apply_mask(alignment, mask):
    return alignment[:, mask]


def _compute_frequencies(alignment):
    return [c.frequencies()
            for c in alignment.iter_positions(ignore_metadata=True)]


def mask(alignment: skbio.TabularMSA, max_gap_frequency: float =0.05,
         min_conservation: float =0.40) -> skbio.TabularMSA:
    # check that parameters are in range
    if max_gap_frequency < 0.0 or max_gap_frequency > 1.0:
        raise ValueError('max_gap_frequency out of range [0.0, 1.0]: %f' %
                         max_gap_frequency)
    if min_conservation < 0.0 or min_conservation > 1.0:
        raise ValueError('min_conservation out of range [0.0, 1.0]: %f' %
                         min_conservation)
    # compute frequencies of all alphabet characters
    frequencies = _compute_frequencies(alignment)
    # compute gap and conservation masks, and then combine them
    sequence_dtype = alignment.dtype
    gap_mask = _compute_gap_mask(frequencies, sequence_dtype,
                                 max_gap_frequency)
    conservation_mask = _compute_conservation_mask(frequencies, sequence_dtype,
                                                   min_conservation)
    combined_mask = gap_mask & conservation_mask
    # apply the mask and return the resulting alignment
    return _apply_mask(alignment, combined_mask)
