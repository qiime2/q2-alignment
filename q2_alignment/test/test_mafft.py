# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio
import unittest

from qiime.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from q2_alignment import mafft


class MafftTests(TestPluginBase):

    package = 'q2_alignment.test'

    def test_mafft(self):
        input_fp = self.get_data_path('unaligned-dna-sequences-1.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')
        exp = skbio.TabularMSA(
            [skbio.DNA('AGGGGGG', metadata={'id': 'seq1', 'description': ''}),
             skbio.DNA('-GGGGGG', metadata={'id': 'seq2', 'description': ''})]
        )
        result = mafft(input_sequences)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

if __name__ == "__main__":
    unittest.main()
