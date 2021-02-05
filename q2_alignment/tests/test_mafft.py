# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import unittest
from unittest.mock import patch, ANY
import subprocess

import skbio
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat, AlignedDNAFASTAFormat
from qiime2.util import redirected_stdio

from q2_alignment import mafft, mafft_add
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
            result = mafft(input_sequences, n_threads='auto')
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

        self.assertIn('a'*250, obs)
        self.assertIn('b'*250, obs)
        self.assertIn('c'*250, obs)

    def test_duplicate_input_ids(self):
        input_fp = self.get_data_path('unaligned-duplicate-ids.fasta')
        input_sequences = DNAFASTAFormat(input_fp, mode='r')

        with self.assertRaisesRegex(ValueError, 'the unaligned.*id1'):
            with redirected_stdio(stderr=os.devnull):
                mafft(input_sequences)

    def test_mafft_parttree_exception(self):
        input_fp = os.path.join(self.temp_dir.name, 'million.fasta')
        with open(input_fp, "w") as f:
            for i in range(0, 1000002):
                f.write('>%d\nAAGCAAGC\n' % i)
        input_sequences = DNAFASTAFormat(input_fp, mode='r')
        with self.assertRaisesRegex(ValueError, '1 million'):
            with redirected_stdio(stderr=os.devnull):
                mafft(input_sequences)


class MafftAddTests(TestPluginBase):
    package = 'q2_alignment.tests'

    def _prepare_sequence_data(self):
        sequences_fp = self.get_data_path('unaligned-dna-sequences-1.fasta')
        sequences = DNAFASTAFormat(sequences_fp, mode='r')
        alignment_fp = self.get_data_path('aligned-dna-sequences-1.fasta')
        alignment = AlignedDNAFASTAFormat(alignment_fp, mode='r')
        exp = skbio.TabularMSA(
            [skbio.DNA('AGGGGG-',
                       metadata={'id': 'aln-seq-1', 'description': ''}),
             skbio.DNA('AGGGGGG',
                       metadata={'id': 'aln-seq-2', 'description': ''}),
             skbio.DNA('AGGGGGG',
                       metadata={'id': 'seq1', 'description': ''}),
             skbio.DNA('-GGGGGG',
                       metadata={'id': 'seq2', 'description': ''})]
        )

        return alignment, sequences, exp

    def test_mafft_add(self):
        alignment, sequences, exp = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_mafft_add_fragments(self):
        alignment, sequences, exp = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences, addfragments=True)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_mafft_add_flags(self):
        alignment, sequences, exp = self._prepare_sequence_data()

        with patch('q2_alignment._mafft.run_command') as patched_run_cmd:
            with patch('q2_alignment._mafft.skbio.TabularMSA.read',
                       return_value=exp):
                _ = mafft_add(alignment, sequences)
                patched_run_cmd.assert_called_with(
                    ["mafft", "--preservecase", "--inputorder", "--thread",
                     "1", "--add", ANY, ANY], ANY)

                _ = mafft_add(alignment, sequences, addfragments=True)
                patched_run_cmd.assert_called_with(
                    ["mafft", "--preservecase", "--inputorder", "--thread",
                     "1", "--addfragments", ANY, ANY], ANY)

    def test_duplicate_input_ids_in_unaligned(self):
        input_fp = self.get_data_path('unaligned-duplicate-ids.fasta')
        sequences = DNAFASTAFormat(input_fp, mode='r')

        alignment, _, _ = self._prepare_sequence_data()

        with self.assertRaisesRegex(ValueError, 'the unaligned.*id1'):
            with redirected_stdio(stderr=os.devnull):
                mafft_add(alignment, sequences)

    def test_duplicate_input_ids_in_aligned(self):
        input_fp = self.get_data_path('aligned-duplicate-ids-1.fasta')
        alignment = DNAFASTAFormat(input_fp, mode='r')

        _, sequences, _ = self._prepare_sequence_data()

        with self.assertRaisesRegex(ValueError, 'the aligned.*id1'):
            with redirected_stdio(stderr=os.devnull):
                mafft_add(alignment, sequences)

    def test_duplicate_input_ids_across_aligned_and_unaligned(self):
        input_fp = self.get_data_path('aligned-duplicate-ids-2.fasta')
        alignment = DNAFASTAFormat(input_fp, mode='r')

        _, sequences, _ = self._prepare_sequence_data()

        with self.assertRaisesRegex(ValueError, 'aligned and unaligned.*seq1'):
            with redirected_stdio(stderr=os.devnull):
                mafft_add(alignment, sequences)

    def test_long_ids_are_not_truncated_unaligned(self):
        input_fp = self.get_data_path('unaligned-long-ids.fasta')
        sequences = DNAFASTAFormat(input_fp, mode='r')

        alignment, _, _ = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences)

        with open(str(result), 'r') as fh:
            obs = fh.read()

        self.assertIn('a'*250, obs)
        self.assertIn('b'*250, obs)
        self.assertIn('c'*250, obs)
        self.assertIn('aln-seq-1', obs)
        self.assertIn('aln-seq-2', obs)

    def test_long_ids_are_not_truncated_aligned(self):
        input_fp = self.get_data_path('aligned-long-ids.fasta')
        alignment = DNAFASTAFormat(input_fp, mode='r')

        _, sequences, _ = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences)

        with open(str(result), 'r') as fh:
            obs = fh.read()

        self.assertIn('a'*250, obs)
        self.assertIn('b'*250, obs)
        self.assertIn('seq1', obs)
        self.assertIn('seq2', obs)


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
