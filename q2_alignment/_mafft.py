# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

from q2_types.feature_data import (DNAFASTAFormat, AlignedDNAFASTAFormat)


def run_command(cmd, output_fp, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
    with open(output_fp, 'w') as output_f:
        subprocess.run(cmd, stdout=output_f, check=True)


def mafft(sequences: DNAFASTAFormat) -> AlignedDNAFASTAFormat:
    result = AlignedDNAFASTAFormat()
    unaligned_fp = str(sequences)
    aligned_fp = str(result)
    cmd = ["mafft", "--preservecase", unaligned_fp]
    run_command(cmd, aligned_fp)
    return result
