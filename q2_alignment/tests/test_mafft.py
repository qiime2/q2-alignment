# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
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

    def _prepare_sequence_data(self):
        input_fp = self.get_data_path('unaligned-dna-sequences-1.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')
        exp = skbio.TabularMSA(
            [skbio.DNA('AGGGGGG', metadata={'id': 'seq1', 'description': ''}),
             skbio.DNA('-GGGGGG', metadata={'id': 'seq2', 'description': ''})]
        )

        return input_sequences, exp

    def test_mafft(self):
        input_sequences, exp = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft(input_sequences)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_multithreaded_mafft(self):
        input_sequences, exp = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft(input_sequences, n_threads=0)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_long_ids_are_not_truncated(self):
        input_fp = self.get_data_path('unaligned-long-ids.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            result = mafft(input_sequences)

        with open(str(result), 'r') as fh:
            obs = fh.read()

        exp_fp = self.get_data_path('aligned-long-ids.fasta')
        with open(exp_fp, 'r') as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)

    def test_duplicate_input_ids(self):
        input_fp = self.get_data_path('unaligned-duplicate-ids.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')

        with self.assertRaisesRegex(ValueError, 'duplicate.*id1'):
            with redirected_stdio(stderr=os.devnull):
                mafft(input_sequences)


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
