# -----------------------------------------------------------------------------
# Copyright (C) Daniel Standage, 2015. It is licensed under the ISC license,
# see LICENSE.txt. Contact: daniel.standage@gmail.com
# -----------------------------------------------------------------------------

"""
Generators for parsing sequence data records in Fasta and Fastq.
"""


def get_parser(fastq=True, paired=True):
    """
    Return a pointer to the correct parsing function based on the specified
    format of the data.
    """
    if fastq:
        if paired:
            return fastq_paired
        else:
            return fastq
    else:
        assert not paired
        return fasta


def check_record(record, fastq=True, paired=True):
    """
    """
    if fastq:
        if paired:
            assert len(record) == 6
        else:
            assert len(record) == 3
    else:
        assert len(record) == 2


def fasta(fp):
    """
    Generator yields sequence records from Fasta files. Stolen shamelessly from
    http://stackoverflow.com/a/7655072/459780.
    """
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, ''.join(seq)


def fastq(fp):
    """
    Generator yields unpaired sequence records from Fastq files. Only supports
    4-line Fastq format.
    """
    linecount = 0
    name, seq, qual = [None] * 3
    for line in fp:
        linecount += 1
        if linecount % 4 == 1:
            name = line.rstrip()
        elif linecount % 4 == 2:
            seq = line.rstrip()
        elif linecount % 4 == 0:
            qual = line.rstrip()
            yield name, seq, qual
            name, seq, qual = [None] * 3


def fastq_paired(fp):
    """
    Generator yields paired sequence records from Fastq files. Only supports
    4-line Fastq format with interleaved pairs.
    """
    linecount = 0
    name1, seq1, qual1, name2, seq2, qual2 = [None] * 6
    for line in fp:
        linecount += 1
        if linecount % 8 == 1:
            name1 = line.rstrip()
        elif linecount % 8 == 2:
            seq1 = line.rstrip()
        elif linecount % 8 == 4:
            qual1 = line.rstrip()
        elif linecount % 8 == 5:
            name2 = line.rstrip()
        elif linecount % 8 == 6:
            seq2 = line.rstrip()
        elif linecount % 8 == 0:
            qual2 = line.rstrip()
            yield name1, seq1, qual1, name2, seq2, qual2
            name1, seq1, qual1, name2, seq2, qual2 = [None] * 6
