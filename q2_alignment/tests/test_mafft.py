# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import unittest
import subprocess

import skbio
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat, AlignedDNAFASTAFormat
from qiime2.util import redirected_stdio

from q2_alignment import mafft
from q2_alignment._mafft import run_command


class MafftTests(TestPluginBase):

    package = 'q2_alignment.tests'

    def test_mafft(self):
        input_fp = self.get_data_path('unaligned-dna-sequences-1.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')
        exp = skbio.TabularMSA(
            [skbio.DNA('AGGGGGG', metadata={'id': 'seq1', 'description': ''}),
             skbio.DNA('-GGGGGG', metadata={'id': 'seq2', 'description': ''})]
        )
        with redirected_stdio(stderr=os.devnull):
            result = mafft(input_sequences)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)


class RunCommandTests(TestPluginBase):

    package = 'q2_alignment.tests'

    def test_failed_run(self):
        input_fp = self.get_data_path('unaligned-dna-sequences-1.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')
        output_alignment = AlignedDNAFASTAFormat()
        unaligned_fp = str(input_sequences)
        aligned_fp = str(output_alignment)
        cmd = ["mafft", "--not-a-real-parameter", unaligned_fp]
        with self.assertRaises(subprocess.CalledProcessError):
            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, aligned_fp)

    def test_failed_run_not_verbose(self):
        input_fp = self.get_data_path('unaligned-dna-sequences-1.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')
        output_alignment = AlignedDNAFASTAFormat()
        unaligned_fp = str(input_sequences)
        aligned_fp = str(output_alignment)
        cmd = ["mafft", "--not-a-real-parameter", unaligned_fp]
        with self.assertRaises(subprocess.CalledProcessError):
            with redirected_stdio(stderr=os.devnull):
                run_command(cmd, aligned_fp, verbose=False)


if __name__ == "__main__":
    unittest.main()
