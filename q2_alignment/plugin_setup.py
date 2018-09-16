# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from textwrap import dedent

from qiime2.plugin import Plugin, Float, Int, Bool, Str, Range, Citations
from q2_types.feature_data import FeatureData, Sequence, AlignedSequence

import q2_alignment

citations = Citations.load('citations.bib', package='q2_alignment')
plugin = Plugin(
    name='alignment',
    version=q2_alignment.__version__,
    website='https://github.com/qiime2/q2-alignment',
    package='q2_alignment',
    description=('This QIIME 2 plugin provides support for generating '
                 'and manipulating sequence alignments.'),
    short_description='Plugin for generating and manipulating alignments.'
)

plugin.methods.register_function(
    function=q2_alignment.mafft,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={'n_threads': Int % Range(0, None),
                'parttree': Bool},
    outputs=[('alignment', FeatureData[AlignedSequence])],
    input_descriptions={'sequences': 'The sequences to be aligned.'},
    parameter_descriptions={
        'n_threads': 'The number of threads. (Use 0 to automatically use all '
                     'available cores)',
        'parttree': 'This flag is required if the number of sequences being '
                    'aligned are larger than 1000000. Disabled by default'},
    output_descriptions={'alignment': 'The aligned sequences.'},
    name='De novo multiple sequence alignment with MAFFT',
    description=("Perform de novo multiple sequence alignment using MAFFT."),
    citations=[citations['katoh2013mafft']]
)

plugin.methods.register_function(
    function=q2_alignment.sina,
    inputs={'sequences': FeatureData[Sequence],
            'reference': FeatureData[AlignedSequence]},
    parameters={'arb_reference': Str,
                'num_references': Int % Range(0, None),
                'kmer_len': Int % Range(5, 25)},
    outputs=[('alignment', FeatureData[AlignedSequence])],
    input_descriptions={
        'sequences': 'The sequences to be aligned.',
        'reference': 'The reference alignment (QZA format).',
    },
    parameter_descriptions={
        'arb_reference': 'The reference alignment (ARB format).',
        'num_references': 'The number of reference sequences used to '
                          'align each sequence.',
        'kmer_len': 'The size of K used in reference search. '
                    'Use 6-8 for sequences <150bp. '
                    'Use 10 for full length 16S. '},
    output_descriptions={'alignment': 'The aligned sequences.'},
    name='Reference based multiple sequence alignment with SINA',
    description=dedent("""
    Perform reference based multiple sequence alignment using SINA.

    Note:
    To use a reference database in ARB format without QZA conversion,
    use the --p-arb-reference parameter instead of --i-reference.
    """),
    citations=[citations['pruesse2012sina']]
)

plugin.methods.register_function(
    function=q2_alignment.mask,
    inputs={'alignment': FeatureData[AlignedSequence]},
    parameters={'max_gap_frequency': Float % Range(0, 1, inclusive_end=True),
                'min_conservation': Float % Range(0, 1, inclusive_end=True)},
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
                 "columns from an alignment. Default min_conservation was "
                 "chosen to reproduce the mask presented in Lane (1991)."),
    citations=[citations['lane1991']]
)
