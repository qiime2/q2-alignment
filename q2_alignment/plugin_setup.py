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
    package='q2_alignment'
)

plugin.methods.register_function(
    function=q2_alignment.mafft,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={},
    outputs=[('alignment', FeatureData[AlignedSequence])],
    name='De novo multiple sequence alignment with MAFFT',
    description=("Perform de novo multiple sequence alignment using MAFFT.")
)

plugin.methods.register_function(
    function=q2_alignment.mask,
    inputs={'alignment': FeatureData[AlignedSequence]},
    parameters={'max_gap_frequency': Float,
                'min_conservation': Float},
    outputs=[('masked_alignment', FeatureData[AlignedSequence])],
    name='Positional conservation and gap filtering.',
    description=("Remove unconserved and highly gapped positions from an "
                 "alignment.")
)
