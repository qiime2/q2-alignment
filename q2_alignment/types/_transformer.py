# ----------------------------------------------------------------------------
# Copyright (c) 2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.util import duplicate
from q2_alignment.plugin_setup import plugin
from q2_types.feature_data import (
    FASTAFormat,
    AlignedProteinSequencesDirectoryFormat,
    AlignedDNASequencesDirectoryFormat,
)


@plugin.register_transformer
def _1(seqs: FASTAFormat) -> AlignedProteinSequencesDirectoryFormat:
    aligned_dir_fmt = AlignedProteinSequencesDirectoryFormat()
    duplicate(
        src=seqs.path,
        dst=aligned_dir_fmt.path / aligned_dir_fmt.file.pathspec
    )
    return aligned_dir_fmt


@plugin.register_transformer
def _2(seqs: FASTAFormat) -> AlignedDNASequencesDirectoryFormat:
    aligned_dir_fmt = AlignedDNASequencesDirectoryFormat()
    duplicate(
        src=seqs.path,
        dst=aligned_dir_fmt.path / aligned_dir_fmt.file.pathspec
    )
    return aligned_dir_fmt
