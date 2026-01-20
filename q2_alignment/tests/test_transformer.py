# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import skbio

from q2_types.feature_data import (
    FASTAFormat,
    AlignedDNASequencesDirectoryFormat,
    AlignedProteinSequencesDirectoryFormat
)

from qiime2.plugin.testing import TestPluginBase


class TestTransformers(TestPluginBase):
    package = 'q2_types.feature_data.tests'

    def test_fasta_format_to_dna_alignment_dir_fmt(self):
        input, obs = self.transform_format(
            FASTAFormat,
            AlignedDNASequencesDirectoryFormat,
            filename='aligned-dna-sequences.fasta'
        )

        exp = skbio.TabularMSA.read(
            str(input),
            format='fasta',
            constructor=skbio.DNA,
        )
        obs = skbio.TabularMSA.read(
            f'{obs}/aligned-dna-sequences.fasta',
            format='fasta',
            constructor=skbio.DNA,
        )

        self.assertEqual(obs, exp)

    def test_fasta_format_to_protein_alignment_dir_fmt(self):
        input, obs = self.transform_format(
            FASTAFormat,
            AlignedProteinSequencesDirectoryFormat,
            filename='aligned-protein-sequences.fasta'
        )

        exp = skbio.TabularMSA.read(
            str(input),
            format='fasta',
            constructor=skbio.Protein,
        )
        obs = skbio.TabularMSA.read(
            f'{obs}/aligned-protein-sequences.fasta',
            format='fasta',
            constructor=skbio.Protein,
        )

        self.assertEqual(obs, exp)
