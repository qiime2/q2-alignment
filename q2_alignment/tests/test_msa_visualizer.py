import os.path

from q2_types.feature_data import AlignedDNAFASTAFormat
from qiime2.plugin.testing import TestPluginBase
from skbio import DNA, TabularMSA

from q2_alignment._msa_visualizer._visualizer import msa_visualizer


class MSAVisualizerTests(TestPluginBase):
    package = "q2_alignment.tests"

    def _prepare_msa(self):
        msa_fp = self.get_data_path("aligned-dna-sequences-1.fasta")
        msa = AlignedDNAFASTAFormat(msa_fp, mode="r")

        exp = TabularMSA(
            [
                DNA("AGGGGGG", metadata={"id": "seq1", "description": ""}),
                DNA("-GGGGGG", metadata={"id": "seq2", "description": ""}),
            ]
        )

        return msa, exp

    def test_msa_visualizer_minimal(self):
        msa, _ = self._prepare_msa()

        with self.temp_dir as output_dir:
            msa_visualizer(output_dir, msa)

            with open(os.path.join(output_dir, "index.html"), "r") as fh:
                observed = fh.read()

            self.assertIn("2", observed)  # expected number of sequences
            self.assertIn("7", observed)  # expected alignment length
            self.assertIn("1/14(7.14%)", observed)  # expected GC content

            # Expected BLAST Urls
            self.assertIn("QUERY=AGGGGG", observed)
            self.assertIn("QUERY=AGGGGGG", observed)
