# -----------------------------------------------------------------------------
# Copyright (C) Daniel Standage, 2015. It is licensed under the ISC license,
# see LICENSE.txt. Contact: daniel.standage@gmail.com
# -----------------------------------------------------------------------------

"""
Utilities for working with complete sequence records from a database of
Fast[aq] sequences. SHA1 hash of each sequence is stored in memory, instead of
the sequence itself. Therefore, these utilities are very space efficient when
working with thousands of long(ish) sequences such as scaffolds and contigs.
The memory savings will not be quite as drastic when working with millions of
short sequences, such as from an RNA-seq experiment.
"""

import sequniq
import hashlib


def uniqseqs(seqdata, trimdefline=False, checkrevcom=False, fastq=True,
             paired=True):
    """
    Given a file of Fast[aq] sequences `seqdata`, retrieve unique sequences.
    Generator function yields complete Fast[aq] records.
    """
    seqs = {}
    parsefunc = sequniq.parse.get_parser(fastq=fastq, paired=paired)
    for record in parsefunc(seqdata):
        sequniq.parse.check_record(record, fastq=fastq, paired=paired)

        seq = record[1]
        if fastq and paired:
            seq += record[4]
        seqsha = hashlib.sha1(seq).hexdigest()

        if seqsha not in seqs:
            if checkrevcom:
                rseqsha = hashlib.sha1(sequniq.revcomp(seq)).hexdigest()
                if rseqsha not in seqs:
                    seqs[seqsha] = 1
                    yield record
            else:
                seqs[seqsha] = 1
                yield record
