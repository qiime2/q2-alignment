# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import types
import tempfile
import subprocess
import os.path

import skbio


def mafft(unaligned_sequences: types.GeneratorType) -> skbio.TabularMSA:
    # TODO use qiime.TemporaryDirectory(), pending merge of qiime2/qiime2#108
    with tempfile.TemporaryDirectory() as tmpdirname:
        unaligned_fp = os.path.join(tmpdirname, 'unaligned.fasta')
        aligned_fp = os.path.join(tmpdirname, 'aligned.fasta')
        # write the unaligned sequences to file
        with open(unaligned_fp, 'w') as unaligned_f:
            for sequence in unaligned_sequences:
                sequence.write(unaligned_f, format='fasta')
        # build the mafft command
        cmd = ["mafft", "--quiet", unaligned_fp]
        # align the sequences and write the output file
        with open(aligned_fp, 'w') as aligned_f:
            subprocess.run(cmd, stdout=aligned_f)
        # read the output file into a TabularMSA and return it
        return skbio.io.read(file=aligned_fp, format='fasta',
                             into=skbio.TabularMSA, constructor=skbio.DNA,
                             lowercase=True)
