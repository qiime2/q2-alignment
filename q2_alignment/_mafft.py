# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

from q2_types.feature_data import (DNAFASTAFormat, AlignedDNAFASTAFormat)


def mafft(sequences: DNAFASTAFormat) -> AlignedDNAFASTAFormat:
    result = AlignedDNAFASTAFormat()
    unaligned_fp = str(sequences)
    aligned_fp = str(result)
    cmd = ["mafft", "--quiet", "--preservecase", unaligned_fp]
    # align the sequences and write the output file
    with open(aligned_fp, 'w') as aligned_f:
        subprocess.run(cmd, stdout=aligned_f)
    return result
