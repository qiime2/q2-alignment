# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
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
from q2_types.feature_data import (
    DNAFASTAFormat,
    AlignedDNAFASTAFormat,
    ProteinFASTAFormat,
    AlignedProteinFASTAFormat,
)
from qiime2.util import redirected_stdio

from q2_alignment import mafft, mafft_add
from q2_alignment._mafft import (
    run_command,
    SequenceType,
    _validate_sequence_pair,
)


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

    def test_mafft_large(self):
        input_sequences, exp = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft(input_sequences, large=True)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)


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

    def _prepare_sequence_data_2(self):
        # for new alignment using `--keeplength` parameter
        sequences_fp = self.get_data_path('unaligned-dna-sequences-2.fasta')
        sequences = DNAFASTAFormat(sequences_fp, mode='r')
        alignment_fp = self.get_data_path('aligned-dna-sequences-2.fasta')
        alignment = AlignedDNAFASTAFormat(alignment_fp, mode='r')
        exp = skbio.TabularMSA(
            [skbio.DNA('AGGG-GGC',
                       metadata={'id': 'aln-seq-1', 'description': ''}),
             skbio.DNA('AGGGTGGC',
                       metadata={'id': 'aln-seq-2', 'description': ''}),
             skbio.DNA('AGGTTGGC',
                       metadata={'id': 'seq-3', 'description': ''}),
             skbio.DNA('AGGATGGC',
                       metadata={'id': 'seq-4', 'description': ''})]
        )

        return alignment, sequences, exp

    def _prepare_sequence_data_3(self):
        # NOT using `--keeplength` parameter. To compare with
        # _prepare_sequence_data_2 output
        sequences_fp = self.get_data_path('unaligned-dna-sequences-2.fasta')
        sequences = DNAFASTAFormat(sequences_fp, mode='r')
        alignment_fp = self.get_data_path('aligned-dna-sequences-2.fasta')
        alignment = AlignedDNAFASTAFormat(alignment_fp, mode='r')
        exp = skbio.TabularMSA(
            [skbio.DNA('AGG--G---GGC',
                       metadata={'id': 'aln-seq-1', 'description': ''}),
             skbio.DNA('AGG--G--TGGC',
                       metadata={'id': 'aln-seq-2', 'description': ''}),
             skbio.DNA('AGG--TTTTGGC',
                       metadata={'id': 'seq-3', 'description': ''}),
             skbio.DNA('AGGTTA--TGGC',
                       metadata={'id': 'seq-4', 'description': ''})]
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

    def test_mafft_add_fragments_keeplength(self):
        alignment, sequences, exp = self._prepare_sequence_data_2()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences, addfragments=True,
                               keeplength=True)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_mafft_add_fragments_no_keeplength(self):
        alignment, sequences, exp = self._prepare_sequence_data_3()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences, addfragments=True,
                               keeplength=False)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_mafft_add_no_keeplength_large(self):
        alignment, sequences, exp = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences, keeplength=False,
                               large=True)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_mafft_add_keeplength_large(self):
        alignment, sequences, exp = self._prepare_sequence_data()

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(alignment, sequences, keeplength=True,
                               large=True)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.DNA)
        self.assertEqual(obs, exp)

    def test_mafft_add_fragments_large(self):
        alignment, sequences, exp = self._prepare_sequence_data()

        with self.assertRaisesRegex(ValueError, '--p-addfragments and.*er.'):
            with redirected_stdio(stderr=os.devnull):
                mafft_add(alignment, sequences, addfragments=True, large=True)

    def test_mafft_add_flags(self):
        alignment, sequences, exp = self._prepare_sequence_data()

        with patch('q2_alignment._mafft.run_command') as patched_run_cmd:
            with patch('q2_alignment._mafft.skbio.TabularMSA.read',
                       return_value=exp):
                _ = mafft_add(alignment, sequences)
                patched_run_cmd.assert_called_with(
                    ["mafft", "--preservecase", "--inputorder", "--thread",
                     "1", "--add", ANY, ANY], ANY, env=None)

                _ = mafft_add(alignment, sequences, addfragments=True)
                patched_run_cmd.assert_called_with(
                    ["mafft", "--preservecase", "--inputorder", "--thread",
                     "1", "--addfragments", ANY, ANY], ANY, env=None)

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

    def test_is_nucleotide_returns_true_for_nucleotide(self):
        assert SequenceType.NUCLEOTIDE.is_nucleotide() is True
        assert SequenceType.NUCLEOTIDE.is_protein() is False

    def test_is_protein_returns_true_for_protein(self):
        assert SequenceType.PROTEIN.is_protein() is True
        assert SequenceType.PROTEIN.is_nucleotide() is False

    def test_validate_sequence_pair_raises(self):
        fake_alignment = AlignedDNAFASTAFormat()
        fake_sequences = ProteinFASTAFormat()

        with self.assertRaisesRegex(TypeError, "Mismatched sequence type"):
            _validate_sequence_pair(fake_alignment, fake_sequences)

    @patch("q2_alignment._mafft._mafft")
    def test_mafft_sets_protein_sequence_type(self, mock_mafft):
        seqs = ProteinFASTAFormat()
        mafft(seqs)
        mock_mafft.assert_called_once()
        args, _ = mock_mafft.call_args
        assert args[-1] == SequenceType.PROTEIN

    @patch("q2_alignment._mafft._mafft")
    def test_mafft_sets_nucleotide_sequence_type(self, mock_mafft):
        seqs = DNAFASTAFormat()
        mafft(seqs)
        mock_mafft.assert_called_once()
        args, _ = mock_mafft.call_args
        assert args[-1] == SequenceType.NUCLEOTIDE

    @patch("q2_alignment._mafft._mafft")
    def test_mafft__add_sets_protein_sequence_type(self, mock_mafft):
        alignment = AlignedProteinFASTAFormat()
        seqs = ProteinFASTAFormat()
        mafft_add(alignment, seqs)
        mock_mafft.assert_called_once()
        args, _ = mock_mafft.call_args
        assert args[-1] == SequenceType.PROTEIN

    @patch("q2_alignment._mafft._mafft")
    def test_mafft_add_sets_nucleotide_sequence_type(self, mock_mafft):
        alignment = AlignedDNAFASTAFormat()
        seqs = DNAFASTAFormat()
        mafft_add(alignment, seqs)
        mock_mafft.assert_called_once()
        args, _ = mock_mafft.call_args
        assert args[-1] == SequenceType.NUCLEOTIDE

    def test_mafft_protein(self):
        input_fp = self.get_data_path('protein-sequences-1.fasta')
        input_sequences = ProteinFASTAFormat(input_fp, mode='r')
        aligned_fp = self.get_data_path('aligned-protein-sequences-1.fasta')
        exp = AlignedProteinFASTAFormat(aligned_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            result = mafft(input_sequences)
        exp = skbio.io.read(str(exp), into=skbio.TabularMSA,
                            constructor=skbio.Protein)
        obs = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.Protein)

        self.assertEqual(obs, exp)

    def test_mafft_add_protein(self):
        sequences_fp = self.get_data_path('protein-sequences-1.fasta')
        input_sequences = ProteinFASTAFormat(sequences_fp, mode='r')
        aligned_fp = self.get_data_path('aligned-protein-sequences-2.fasta')
        input_alignment = AlignedProteinFASTAFormat(aligned_fp, mode='r')
        exp_fp = self.get_data_path('aligned-protein-sequences-3.fasta')
        exp_alignment = AlignedProteinFASTAFormat(exp_fp, mode='r')

        with redirected_stdio(stderr=os.devnull):
            result = mafft_add(input_alignment, input_sequences)
        exp = skbio.io.read(str(result), into=skbio.TabularMSA,
                            constructor=skbio.Protein)
        obs = skbio.io.read(str(exp_alignment), into=skbio.TabularMSA,
                            constructor=skbio.Protein)

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
