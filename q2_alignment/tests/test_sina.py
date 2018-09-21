# ----------------------------------------------------------------------------
# Copyright (c) 2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

__author__ = "Elmar Pruesse"

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
        """Prepares test data tuples: reference, query and expected

        This function generates 1) reference alignments as subsets of length
        n-1 from the reference alignment found in "aligned_fna.fasta", a
        query (the removed sequence w/o gaps) and an expected alignment
        (the removed sequence w gaps).
        """
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
        """Tests generally correct function of sina() method

        We expect that when re-aligning each individual sequence in our test
        data set of 18 short 16S fragments (using the remaining 17 sequences as
        a reference), the original alignment is reproduced with at most two
        errors, which amounts to >99.5% identity with the original alignment.
        (It's actually 100% right now, but this test should not be fragile
        w.r.t SINA versions).
        """
        count = 0
        sum_match_f = 0
        for ref, query, exp in self._prepare_sequence_data():
            with redirected_stdio(stderr=os.devnull):
                result = sina(query, ref, kmer_len=6, num_references=5)
            aligned = TabularMSA.read(str(result), constructor=DNA)
            sum_match_f += aligned[0].match_frequency(exp, relative=True)
            count += 1
        avg_match_frequency = sum_match_f / count
        self.assertTrue(avg_match_frequency > 0.995)

    def test_duplicate_input_ids(self):
        input_fp = self.get_data_path('aligned_dna.fasta')
        input_sequences = AlignedDNAFASTAFormat(input_fp, mode='r')
        ref = AlignedDNAFASTAFormat()
        msa = TabularMSA.read(str(input_sequences), constructor=DNA)
        ref_msa = msa[0:3]
        ref_msa.extend(msa[0:1], index=msa.index[4:5])
        ref_msa.write(ref.open())
        query = DNAFASTAFormat()
        msa[6].degap().write(query.open())

        with self.assertRaisesRegex(ValueError, 'Duplicate.*AceElong'):
            with redirected_stdio(stderr=os.devnull):
                sina(query, ref)

    def test_params(self):
        ref, query, exp = next(self._prepare_sequence_data())
        with self.assertRaisesRegex(ValueError, 'Only either'):
            with redirected_stdio(stderr=os.devnull):
                sina(query, ref, str(ref))
        with self.assertRaisesRegex(ValueError, 'needs a reference'):
            with redirected_stdio(stderr=os.devnull):
                sina(query)

if __name__ == "__main__":
    unittest.main()
