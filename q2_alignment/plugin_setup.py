# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin, Float
from q2_types.feature_data import FeatureData, Sequence, AlignedSequence

import q2_alignment

plugin = Plugin(
    name='alignment',
    version=q2_alignment.__version__,
    website='https://github.com/qiime2/q2-alignment',
    package='q2_alignment',
    citation_text=("MAFFT multiple sequence alignment software version 7: "
                   "improvements in performance and usability. Katoh K, "
                   "Standley DM. Mol Biol Evol. 2013 Apr;30(4):772-80. "
                   "doi: 10.1093/molbev/mst010."),
    description=('This QIIME 2 plugin provides support for generating '
                 'and working with sequence alignments.'),
    short_description='Plugin for generating and manipulating alignments.'
)

plugin.methods.register_function(
    function=q2_alignment.mafft,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={},
    outputs=[('alignment', FeatureData[AlignedSequence])],
    input_descriptions={'sequences': 'The sequences to be aligned.'},
    parameter_descriptions={},
    output_descriptions={'alignment': 'The aligned sequences.'},
    name='De novo multiple sequence alignment with MAFFT',
    description=("Perform de novo multiple sequence alignment using MAFFT.")
)

plugin.methods.register_function(
    function=q2_alignment.mask,
    inputs={'alignment': FeatureData[AlignedSequence]},
    parameters={'max_gap_frequency': Float,
                'min_conservation': Float},
    outputs=[('masked_alignment', FeatureData[AlignedSequence])],
    input_descriptions={'alignment': 'The alignment to be masked.'},
    parameter_descriptions={
        'max_gap_frequency': ('The maximum relative frequency of gap '
                              'characters in a column for the column to be '
                              'retained. This relative frequency must be a '
                              'number between 0.0 and 1.0 (inclusive), where '
                              '0.0 retains only those columns without gap '
                              'characters, and 1.0 retains all columns '
                              'regardless of gap character frequency.'),
        'min_conservation': ('The minimum relative frequency '
                             'of at least one non-gap character in a '
                             'column for that column to be retained. This '
                             'relative frequency must be a number between 0.0 '
                             'and 1.0 (inclusive). For example, if a value of '
                             '0.4 is provided, a column will only be retained '
                             'if it contains at least one character that is '
                             'present in at least 40% of the sequences.')
    },
    output_descriptions={'masked_alignment': 'The masked alignment.'},
    name='Positional conservation and gap filtering.',
    description=("Mask (i.e., filter) unconserved and highly gapped "
                 "columns from an alignment. ")
)
