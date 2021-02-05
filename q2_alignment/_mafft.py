# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

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


def _mafft(sequences_fp, alignment_fp, n_threads, parttree, addfragments):
    # Save original sequence IDs since long ids (~250 chars) can be truncated
    # by mafft. We'll replace the IDs in the aligned sequences file output by
    # mafft with the originals.
    #
    # https://github.com/qiime2/q2-alignment/issues/37
    aligned_seq_ids = {}
    unaligned_seq_ids = {}

    if alignment_fp is not None:
        for seq in skbio.io.read(alignment_fp, format='fasta',
                                 constructor=skbio.DNA):
            id_ = seq.metadata['id']
            if id_ in aligned_seq_ids:
                raise ValueError(
                    "A sequence ID is duplicated in the aligned sequences: "
                    "%r" % id_)
            else:
                aligned_seq_ids[id_] = True

    for seq in skbio.io.read(sequences_fp, format='fasta',
                             constructor=skbio.DNA):
        id_ = seq.metadata['id']
        if id_ in unaligned_seq_ids:
            raise ValueError(
                "A sequence ID is duplicated in the unaligned sequences: "
                "%r" % id_)
        elif id_ in aligned_seq_ids:
            raise ValueError(
                "A sequence ID is present in both the aligned and unaligned "
                "sequences: %r" % id_)
        else:
            unaligned_seq_ids[id_] = True

    result = AlignedDNAFASTAFormat()
    result_fp = str(result)
    ids = {**aligned_seq_ids, **unaligned_seq_ids}

    # mafft will fail if the number of sequences is larger than 1 million.
    # mafft requires using parttree which is an algorithm to build an
    # approximate tree from a large number of unaligned sequences.
    # By catching the error below if a user has not used parttree flag, we are
    # eliminating the need for the mafft error to be shown to the user which
    # can be confusing and intimidating.

    if not parttree and len(ids) > 1000000:
        raise ValueError(
            "The number of sequences in your feature table is larger than "
            "1 million, please use the parttree parameter")

    # mafft's signal for utilizing all cores is -1. We want to our users
    # to enter auto for using all cores. This is to prevent any confusion and
    # to keep the UX consisent.
    if n_threads == 'auto':
        n_threads = -1

    # `--inputorder` must be turned on because we need the input and output in
    # the same sequence order to replace the IDs below. This is mafft's default
    # behavior but we pass the flag in case that changes in the future.
    cmd = ["mafft", "--preservecase", "--inputorder",
           "--thread", str(n_threads)]

    if parttree:
        cmd += ['--parttree']

    if alignment_fp is not None:
        add_flag = '--addfragments' if addfragments else '--add'
        cmd += [add_flag, sequences_fp, alignment_fp]
    else:
        cmd += [sequences_fp]

    run_command(cmd, result_fp)

    # Read output alignment into memory, reassign original sequence IDs, and
    # write alignment back to disk.
    msa = skbio.TabularMSA.read(result_fp, format='fasta',
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
    msa.write(result_fp, id_whitespace_replacement=None,
              description_newline_replacement=None)
    return result


def mafft(sequences: DNAFASTAFormat,
          n_threads: int = 1,
          parttree: bool = False) -> AlignedDNAFASTAFormat:
    sequences_fp = str(sequences)
    return _mafft(sequences_fp, None, n_threads, parttree, False)


def mafft_add(alignment: AlignedDNAFASTAFormat,
              sequences: DNAFASTAFormat,
              n_threads: int = 1,
              parttree: bool = False,
              addfragments: bool = False) -> AlignedDNAFASTAFormat:
    alignment_fp = str(alignment)
    sequences_fp = str(sequences)
    return _mafft(
        sequences_fp, alignment_fp, n_threads, parttree, addfragments)
