# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import skbio
import unittest

from q2_alignment import mafft


class MafftTests(unittest.TestCase):

    def test_mafft(self):
        input_sequences = [skbio.DNA('AGGGGGG', metadata={'id': 'seq1'}),
                           skbio.DNA('GGGGGG', metadata={'id': 'seq2'})]
        output_alignment = skbio.TabularMSA(
            [skbio.DNA('AGGGGGG', metadata={'id': 'seq1', 'description': ''}),
             skbio.DNA('-GGGGGG', metadata={'id': 'seq2', 'description': ''})]
        )
        result = mafft(input_sequences)
        self.assertEqual(result, output_alignment)

if __name__ == "__main__":
    unittest.main()
