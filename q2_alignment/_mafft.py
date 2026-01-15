# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from enum import Enum
import os
import subprocess
from typing import Union

import skbio
import skbio.io
from q2_types.feature_data import (
    DNAFASTAFormat,
    AlignedDNAFASTAFormat,
    ProteinFASTAFormat,
    AlignedProteinFASTAFormat,
    FASTAFormat
)
from qiime2 import get_cache


class SequenceType(str, Enum):
    NUCLEOTIDE = "nucleotide"
    PROTEIN = "protein"

    def is_nucleotide(self) -> bool:
        return self is SequenceType.NUCLEOTIDE

    def is_protein(self) -> bool:
        return self is SequenceType.PROTEIN


def _validate_sequence_pair(alignment, sequences):
    if (
        isinstance(alignment, AlignedDNAFASTAFormat)
        and not isinstance(sequences, DNAFASTAFormat)
    ) or (
        isinstance(alignment, AlignedProteinFASTAFormat)
        and not isinstance(sequences, ProteinFASTAFormat)
    ):
        raise TypeError(
            "Mismatched sequence types: 'alignment' and 'sequences' must both "
            "be either DNA or protein."
        )


def run_command(cmd, output_fp, verbose=True, env=None):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
    with open(output_fp, 'w') as output_f:
        subprocess.run(cmd, stdout=output_f, check=True, env=env)


def _mafft(sequences_fp, alignment_fp, n_threads, parttree, addfragments,
           keeplength, large, strategy, maxiterate, retree, sequence_type):
    # Save original sequence IDs since long ids (~250 chars) can be truncated
    # by mafft. We'll replace the IDs in the aligned sequences file output by
    # mafft with the originals.
    #
    # https://github.com/qiime2/q2-alignment/issues/37
    aligned_seq_ids = {}
    unaligned_seq_ids = {}

    constructor = skbio.DNA if sequence_type.is_nucleotide() else skbio.Protein

    if alignment_fp is not None:
        for seq in skbio.io.read(alignment_fp, format='fasta',
                                 constructor=constructor):
            id_ = seq.metadata['id']
            if id_ in aligned_seq_ids:
                raise ValueError(
                    "A sequence ID is duplicated in the aligned sequences: "
                    "%r" % id_)
            else:
                aligned_seq_ids[id_] = True

    for seq in skbio.io.read(sequences_fp, format='fasta',
                             constructor=constructor):
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

    result = FASTAFormat()
    result.aligned = True
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

    env = None

    # mafft's signal for utilizing all cores is -1. We want to our users
    # to enter auto for using all cores. This is to prevent any confusion and
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

    if keeplength:
        cmd += ['--keeplength']

    if large and addfragments:
        raise ValueError('--p-addfragments and --p-large cannot be used '
                         'together.')
    elif large:
        env = os.environ.copy()
        env.update({'MAFFT_TMPDIR': get_cache().get_tmp_path()})
        cmd += ['--large']

    if strategy:
        cmd += [("--" + strategy)]

    if maxiterate is not None:
        cmd += ['--maxiterate', str(maxiterate)]

    if retree is not None:
        cmd += ['--retree', str(retree)]

    if alignment_fp is not None:
        add_flag = '--addfragments' if addfragments else '--add'
        cmd += [add_flag, sequences_fp, alignment_fp]

    else:
        cmd += [sequences_fp]

    run_command(cmd, result_fp, env=env)

    # Read output alignment into memory, reassign original sequence IDs, and
    # write alignment back to disk.
    msa = skbio.TabularMSA.read(result_fp, format='fasta',
                                constructor=constructor)
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


def mafft(sequences: Union[DNAFASTAFormat, ProteinFASTAFormat],
          n_threads: int = 1,
          parttree: bool = False,
          large: bool = False,
          strategy: str | None = None,
          maxiterate: int | None = None,
          retree: int | None = None,) -> FASTAFormat:
    sequence_type = SequenceType.NUCLEOTIDE
    if isinstance(sequences, ProteinFASTAFormat):
        sequence_type = SequenceType.PROTEIN

    sequences_fp = str(sequences)

    return _mafft(
        sequences_fp, None, n_threads, parttree, False, False, large,
        strategy, maxiterate, retree, sequence_type)


def mafft_add(alignment: Union[AlignedDNAFASTAFormat,
                               AlignedProteinFASTAFormat],
              sequences: Union[DNAFASTAFormat, ProteinFASTAFormat],
              n_threads: int = 1,
              parttree: bool = False,
              addfragments: bool = False,
              keeplength: bool = False,
              large: bool = False,
              strategy: str | None = None,
              maxiterate: int | None = None,
              retree: int | None = None) -> FASTAFormat:
    _validate_sequence_pair(alignment, sequences)

    sequence_type = SequenceType.NUCLEOTIDE
    if isinstance(sequences, ProteinFASTAFormat):
        sequence_type = SequenceType.PROTEIN

    alignment_fp = str(alignment)
    sequences_fp = str(sequences)

    return _mafft(
        sequences_fp, alignment_fp, n_threads, parttree, addfragments,
        keeplength, large, strategy, maxiterate, retree, sequence_type)
