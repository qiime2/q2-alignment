# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.plugin import Plugin
from q2_types import FeatureData, Sequence, AlignedSequence

import q2_alignment

plugin = Plugin(
    name='alignment',
    version=q2_alignment.__version__,
    website='https://github.com/qiime2/q2-alignment',
    package='q2_alignment'
)

plugin.methods.register_function(
    function=q2_alignment.mafft,
    inputs={'unaligned_sequences': FeatureData[Sequence]},
    parameters={},
    outputs=[('aligned_sequences', FeatureData[AlignedSequence])],
    name='De novo multiple sequence alignment with MAFFT',
    description=("Perform de novo multiple sequence alignment using MAFFT.")
)
