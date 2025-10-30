# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Plugin, Float, Bool, Range, Citations, Threads, Int)
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
    parameters={'n_threads': Threads,
                'parttree': Bool,
                'large': Bool,
                'globalpair': Bool,
                'localpair': Bool,
                'genafpair': Bool,
                'maxiterate': Int,
                'retree': Int,
                'nofft': Bool,
                'auto': Bool},
    outputs=[('alignment', FeatureData[AlignedSequence])],
    input_descriptions={'sequences': 'The sequences to be aligned.'},
    parameter_descriptions={
        'n_threads': 'The number of threads. (Use `auto` to automatically use '
                     'all available cores)',
        'parttree': 'This flag is required if the number of sequences being '
                    'aligned are larger than 1000000. Disabled by default.',
        'large': 'This flag is required when aligning very large datasets '
                 'that do not otherwise fit into memory. Temporary data is '
                 'then stored in files, instead of RAM. The --use-cache '
                 'flag specifies the storage location of the temporary files '
                 'created. By default, $TMP/qiime2/ is used.',
        'globalpair': 'Compute all pairwise alignments using the '
                      'Needleman-Wunsch algorithm. Suitable for up to ~200 '
                      'sequences. A combination with --p-maxiterate 1000 '
                      'is recommended (G-INS-i).',
        'localpair': 'Compute all pairwise alignments using the '
                     'Smith-Waterman algorithm. Suitable for up to ~200 '
                     'sequences. A combination with --p-maxiterate 1000 '
                     'is recommended (L-INS-i).',
        'genafpair': 'Compute all pairwise alignments with a local algorithm '
                     'with the generalized affine gap cost. Suitable for up '
                     'to ~200 sequences. A combination with --p-maxiterate '
                     '1000 is recommended (E-INS-i).',
        'maxiterate': 'Specifies how many iterative refinement cycles are '
                      'performed after the initial progressive alignment. '
                      'By default, no iterative refinement is performed.',
        'retree': 'Specifies the number of times the guide tree is rebuilt '
                  'during the progressive stage. Typically, tree topology '
                  'stabilizes after 2-3 iterations and higher values rarely '
                  'improves alignment quality enough to justify the extra '
                  'computation.',
        'nofft': 'Disables Fast Fourier Transform (FFT) approximation in '
                 'group-to-group alignment. In general, the FFT algorithm is '
                 'less efficient for distantly related (i.e., less '
                 'conserved) sequences.',
        'auto': 'Automatically select the best alignment strategy '
                '(from FFT-NS-1, FFT-NS-2, FFT-NS-i, or L-INS-i) based on '
                'the input\'s data size.'},
    output_descriptions={'alignment': 'The aligned sequences.'},
    name='De novo multiple sequence alignment with MAFFT',
    description=("Perform de novo multiple sequence alignment using MAFFT."),
    citations=[citations['katoh2013mafft']]
)

plugin.methods.register_function(
    function=q2_alignment.mafft_add,
    inputs={'alignment': FeatureData[AlignedSequence],
            'sequences': FeatureData[Sequence]},
    parameters={'n_threads': Threads,
                'parttree': Bool,
                'addfragments': Bool,
                'keeplength': Bool,
                'large': Bool,
                'globalpair': Bool,
                'localpair': Bool,
                'genafpair': Bool,
                'maxiterate': Int,
                'retree': Int,
                'nofft': Bool,
                'auto': Bool},
    outputs=[('expanded_alignment', FeatureData[AlignedSequence])],
    input_descriptions={'alignment': 'The alignment to which '
                                     'sequences should be added.',
                        'sequences': 'The sequences to be added.'},
    parameter_descriptions={
        'n_threads': 'The number of threads. (Use `auto` to automatically use '
                     'all available cores)',
        'parttree': 'This flag is required if the number of sequences being '
                    'aligned are larger than 1000000. Disabled by default.',
        'addfragments': 'Optimize for the addition of short sequence '
                        'fragments (for example, primer or amplicon '
                        'sequences). If not set, default sequence addition '
                        'is used.',
        'keeplength': 'If selected, the alignment length will be unchanged. '
                      'Any added sequence that would otherwise introduce new '
                      'insertions into the alignment, will have those '
                      'insertions deleted, to preserve original alignment '
                      'length.',
        'large': 'This flag is required when aligning very large datasets '
                 'that do not otherwise fit into memory. Temporary data is '
                 'then stored in files, instead of RAM. The --use-cache '
                 'flag specifies the storage location of the temporary files '
                 'created. By default, $TMP/qiime2/ is used.',
        'globalpair': 'Compute all pairwise alignments using the '
                      'Needleman-Wunsch algorithm. Suitable for up to ~200 '
                      'sequences. A combination with --p-maxiterate 1000 '
                      'is recommended (G-INS-i).',
        'localpair': 'Compute all pairwise alignments using the Smith-Waterman'
                     'algorithm. Suitable for up to ~200 sequences. A '
                     'combination with --p-maxiterate 1000 is recommended '
                     '(L-INS-i).',
        'genafpair': 'Compute all pairwise alignments with a local algorithm '
                     'with the generalized affine gap cost. Suitable for up '
                     'to ~200 sequences. A combination with --p-maxiterate '
                     '1000 is recommended (E-INS-i).',
        'maxiterate': 'Specifies how many iterative refinement cycles are '
                      'performed after the initial progressive alignment. '
                      'By default, no iterative refinement is performed.',
        'retree': 'Specifies the number of times the guide tree is rebuilt '
                  'during the progressive stage. Typically, tree topology '
                  'stabilizes after 2-3 iterations and higher values rarely '
                  'improves alignment quality enough to justify the extra '
                  'computation.',
        'nofft': 'Disables Fast Fourier Transform (FFT) approximation in '
                 'group-to-group alignment. In general, the FFT algorithm is '
                 'less efficient for distantly related (i.e., less '
                 'conserved) sequences.',
        'auto': 'Automatically select the best alignment strategy '
                '(from FFT-NS-1, FFT-NS-2, FFT-NS-i, or L-INS-i) based on '
                'the input\'s data size.'},
    output_descriptions={
        'expanded_alignment': 'Alignment containing the provided aligned and '
                              'unaligned sequences.'},
    name='Add sequences to multiple sequence alignment with MAFFT.',
    description='Add new sequences to an existing alignment with MAFFT.',
    citations=[citations['katoh2013mafft']]
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
