# ----------------------------------------------------------------------------
# Copyright (c) 2018, Elmar Pruesse.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path as op
import subprocess as sp
from tempfile import TemporaryDirectory

import click

from q2_types.feature_data import AlignedDNAFASTAFormat, DNAFASTAFormat
from qiime2.plugin import Int, Str

import skbio
import skbio.io


def _run_command(cmd):
    print("# Running command: " + " ".join(cmd))
    sp.run(cmd, check=True)
    print("# Command "+cmd[0]+" finished")


def check_no_duplicate_fasta_id(fasta_fp):
    ids = set()
    for seq in skbio.io.read(fasta_fp, format='fasta',
                             constructor=skbio.DNA):
        fasta_id = seq.metadata['id']
        if fasta_id in ids:
            raise ValueError(
                "Duplicate FastA ID '%' in file '%'"
                % (id))
        ids.add(fasta_id)


def sina(sequences: DNAFASTAFormat,
         reference: AlignedDNAFASTAFormat=None,
         arb_reference: Str=None,
         num_references: Int=40,
         kmer_len: Int=10) -> AlignedDNAFASTAFormat:
    if not reference and not arb_reference:
        raise ValueError(
            "SINA needs a reference alignment.\n\n"
            "Please use either\n"
            "    '--i-reference=reference.qza'\n"
            "  or\n"
            "    '--p-arb-reference=reference.arb'\n"
            "to indicate the reference alignent you wish to use."
        )
    if reference and arb_reference:
        raise ValueError(
            "Only either -i-reference or --p-arb-reference may be specified"
        )

    aligned = AlignedDNAFASTAFormat()

    with TemporaryDirectory() as tmpdir:
        if not arb_reference:  # Convert QZA aligned FAST to ARB
            check_no_duplicate_fasta_id(fasta_fp)
            arb_reference = op.join(tmpdir, "reference.arb")
            _run_command([
                "sina",
                "--intype", "FASTA",
                "--in", str(reference),
                "--outtype", "ARB",
                "--out", arb_reference,
                "--prealigned",
            ])

        _run_command([
            "sina",
            "--intype", "FASTA",
            "--in", str(sequences),
            "--outtype", "FASTA",
            "--out", str(aligned),
            "--fasta-write-dna",
            "--db", arb_reference,
            "--fs-req-gaps", "1",  # q2 checks ref is aligned
            "--fs-min-len", "10",  # allow short references
            "--fs-kmer-len", str(kmer_len),  # set k
            "--fs-max", str(num_references),
        ])

    return aligned
