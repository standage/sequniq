# -----------------------------------------------------------------------------
# Copyright (C) Daniel Standage, 2015. It is licensed under the ISC license,
# see LICENSE.txt. Contact: daniel.standage@gmail.com
# -----------------------------------------------------------------------------

"""
Utilities for working with IDs of unique and non-unique sequences in a database
of Fast[aq] sequences. SHA1 hash of each sequence is stored in memory, instead
of the sequence itself. Therefore, these utilities are very space efficient
when working with thousands of long(ish) sequences such as scaffolds and
contigs. The memory savings will not be quite as drastic when working with
millions of short sequences, such as from an RNA-seq experiment.
"""

import sequniq
import hashlib


def group_seqids_by_sha1(seqdata, trimdefline=True, checkrevcom=True,
                         fastq=False, paired=False):
    """
    Given a file of Fast[aq] sequences `seqdata`, create a dictionary that
    associates each sequence with its SHA1 hash (or the hash of its reverse
    complement). Dictionary key is the hash; dictionary value is a list with 1
    or more sequence IDs.
    """
    seqids = {}
    parsefunc = sequniq.parse.get_parser(fastq=fastq, paired=paired)
    for record in parsefunc(seqdata):
        sequniq.parse.check_record(record, fastq=fastq, paired=paired)

        seqid = record[0]
        if trimdefline:
            seqid = seqid[1:].split(" ")[0]
        seq = record[1]
        if fastq and paired:
            seq += record[4]

        seqsha = hashlib.sha1(seq).hexdigest()
        if seqsha in seqids:
            seqids[seqsha].append(seqid)
        else:
            rseqsha = None
            if checkrevcom:
                rseqsha = hashlib.sha1(sequniq.revcomp(seq)).hexdigest()
            if rseqsha and rseqsha in seqids:
                seqids[rseqsha].append(seqid)
            else:
                seqids[seqsha] = []
                seqids[seqsha].append(seqid)

    return seqids


def dupids(seqdata, trimdefline=True, checkrevcom=True, fastq=False,
           paired=False):
    """
    Given a file of Fast[aq] sequences `seqdata`, determine duplicate sequences
    and provide their IDs. Generator function yields lists of sequence IDs.
    """
    seqids = group_seqids_by_sha1(
        seqdata,
        trimdefline=trimdefline,
        checkrevcom=checkrevcom,
        fastq=fastq,
        paired=paired)

    for seqsha in seqids:
        idlist = seqids[seqsha]
        if len(idlist) > 1:
            yield idlist


def uniqids(seqdata, trimdefline=True, checkrevcom=True, fastq=False,
            paired=False):
    """
    Given a file of Fasta sequences `seqdata`, determine unique sequences and
    provide their IDs. Generator function yields lists of sequence IDs.
    """
    seqids = group_seqids_by_sha1(
        seqdata,
        trimdefline=trimdefline,
        checkrevcom=checkrevcom,
        fastq=fastq,
        paired=paired)

    for seqsha in seqids:
        idlist = seqids[seqsha]
        if len(idlist) > 1:
            yield idlist
