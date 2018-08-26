# ----------------------------------------------------------------------------
# Copyright (c) 2018, Elmar Pruesse.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess as sp

from q2_types.feature_data import AlignedDNAFASTAFormat, DNAFASTAFormat

import skbio
import skbio.io


def sina(sequences: DNAFASTAFormat,
         reference: DNAFASTAFormat,
         n_threads: int=1) -> AlignedDNAFASTAFormat:
    unaligned_fp = str(sequences)
    reference_fp = str(reference)
    aligned = AlignedDNAFASTAFormat()
    aligned_fp = str(aligned)

    # Guard against duplicate IDs
    ids = set()
    for seq in skbio.io.read(unaligned_fp, format='fasta',
                             constructor=skbio.DNA):
        fasta_id = seq.metadata['id']
        if fasta_id in ids:
            raise ValueError(
                "Encountered duplicate sequence ID in unaligned sequences: %r"
                % id)
        ids.add(fasta_id)

    cmd = [
        "sina",
        "--intype", "FASTA"
        "--in", unaligned_fp,
        "--outtype", "FASTA",
        "--out", aligned_fp,
        "--ptdb", reference_fp
    ]

    sp.run(cmd, check=True)

    return aligned
