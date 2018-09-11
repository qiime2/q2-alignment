# ----------------------------------------------------------------------------
# Copyright (c) 2018, Elmar Pruesse.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import unittest

from q2_alignment import sina

from q2_types.feature_data import AlignedDNAFASTAFormat, DNAFASTAFormat

from qiime2.plugin.testing import TestPluginBase
from qiime2.util import redirected_stdio

from skbio import TabularMSA, DNA


class SINATests(TestPluginBase):
    package = 'q2_alignment.tests'

    def _prepare_sequence_data(self):
        input_fp = self.get_data_path('aligned_dna.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        msa = TabularMSA.read(str(input_sequences), constructor=DNA)
        for n in range(len(msa)):
            ref = AlignedDNAFASTAFormat()
            ref_msa = msa[0:n]
            ref_msa.extend(msa[n+1:len(msa)], index=msa.index[n+1:len(msa)])
            ref_msa.write(ref.open())
            query = DNAFASTAFormat()
            msa[n].degap().write(query.open())
            exp = msa[n]
            yield ref, query, exp

    def test_sina(self):
        for ref, query, exp in self._prepare_sequence_data():
            with redirected_stdio(stderr=os.devnull):
                result = sina(query, ref, kmer_len=6)
            aligned = TabularMSA.read(str(result), constructor=DNA)
            self.assertEqual(aligned[0], exp)


if __name__ == "__main__":
    unittest.main()
