# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import subprocess

import skbio
import skbio.io
from q2_types.feature_data import DNAFASTAFormat, AlignedDNAFASTAFormat


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


def mafft(sequences: DNAFASTAFormat,
          n_threads: int=1,
          parttree: bool=False) -> AlignedDNAFASTAFormat:
    unaligned_fp = str(sequences)

    # Save original sequence IDs since long ids (~250 chars) can be truncated
    # by mafft. We'll replace the IDs in the aligned sequences file output by
    # mafft with the originals.
    #
    # https://github.com/qiime2/q2-alignment/issues/37
    #
    # Note: using OrderedDict to maintain order of IDs and have quick lookup
    # for duplicates.
    ids = collections.OrderedDict()
    for seq in skbio.io.read(unaligned_fp, format='fasta',
                             constructor=skbio.DNA):
        id = seq.metadata['id']
        if id in ids:
            raise ValueError(
                "Encountered duplicate sequence ID in unaligned sequences: %r"
                % id)
        else:
            ids[id] = True

    result = AlignedDNAFASTAFormat()
    aligned_fp = str(result)

    if not parttree and len(ids) > 1000000:
        raise ValueError(
            "The number of sequences in your feature table is larger than "
            "1 million, please use the parttree parameter")

    # mafft's signal for utilizing all cores is -1. We want to our users
    # to enter 0 for using all cores. This is to prevent any confusion and
    # to keep the UX consisent.
    if n_threads == 0:
        n_threads = -1

    # `--inputorder` must be turned on because we need the input and output in
    # the same sequence order to replace the IDs below. This is mafft's default
    # behavior but we pass the flag in case that changes in the future.
    cmd = ["mafft", "--preservecase", "--inputorder",
           "--thread", str(n_threads)]

    if parttree:
        cmd += ['--parttree']
    
    cmd += ['unaligned_fp']
    run_command(cmd, aligned_fp)

    # Read output alignment into memory, reassign original sequence IDs, and
    # write alignment back to disk.
    msa = skbio.TabularMSA.read(aligned_fp, format='fasta',
                                constructor=skbio.DNA)
    # Using `assert` because mafft would have had to add or drop sequences
    # while aligning, which would be a bug on mafft's end. This is just a
    # sanity check and is not expected to trigger in practice.
    assert len(ids) == len(msa)
    for id, seq in zip(ids, msa):
        seq.metadata['id'] = id

    # Turning off roundtripping options to speed up writing. We can safely turn
    # these options off because we know the sequence IDs are rountrip-safe
    # since we read them from a FASTA file above.
    #
    # http://scikit-bio.org/docs/latest/generated/
    #     skbio.io.format.fasta.html#writer-specific-parameters
    msa.write(aligned_fp, id_whitespace_replacement=None,
              description_newline_replacement=None)
    return result
