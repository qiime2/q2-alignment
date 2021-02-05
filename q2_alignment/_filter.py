# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
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
    # this try/accept allows us to support the ignore_metadata optimization
    # available in scikit-bio 0.5.0 and later, but be backward compatible
    # with scikit-bio 0.4.2, which some plugins depend on.
    try:
        return [c.frequencies()
                for c in alignment.iter_positions(ignore_metadata=True)]
    except TypeError:
        return [c.frequencies()
                for c in alignment.iter_positions()]


def mask(alignment: skbio.TabularMSA, max_gap_frequency: float = 1.0,
         min_conservation: float = 0.40) -> skbio.TabularMSA:
    # check that parameters are in range
    if max_gap_frequency < 0.0 or max_gap_frequency > 1.0:
        raise ValueError('max_gap_frequency out of range [0.0, 1.0]: %f' %
                         max_gap_frequency)
    if min_conservation < 0.0 or min_conservation > 1.0:
        raise ValueError('min_conservation out of range [0.0, 1.0]: %f' %
                         min_conservation)
    # check that input alignment is not empty
    if alignment.shape.position == 0:
        raise ValueError('Input alignment is empty (i.e., there are zero '
                         'sequences or positions in the input alignment).')
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
    result = _apply_mask(alignment, combined_mask)

    if result.shape.position == 0:
        num_input_positions = alignment.shape.position
        frac_passed_gap = (gap_mask.sum() / num_input_positions)
        str_passed_gap = '{percent:.2%}'.format(percent=frac_passed_gap)
        frac_passed_conservation = \
            (conservation_mask.sum() / num_input_positions)
        str_passed_conservation = \
            '{percent:.2%}'.format(percent=frac_passed_conservation)
        raise ValueError("No alignment positions remain after filtering. The "
                         "filter thresholds will need to be relaxed. %s "
                         "of positions were retained by the gap filter, and "
                         "%s of positions were retained by the "
                         "conservation filter." %
                         (str_passed_gap, str_passed_conservation))
    return result
